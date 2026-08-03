# R/ms1.R — MS1 EIC → peak → metric pipeline (+ EIC/peaks/stats/plots helpers).
#
# Pipeline stages (all backend-agnostic except load_spectra_ms1):
#   1. load_spectra_ms1(ds, files) → MS1 Spectra
#   2. derive_extraction_config(...) → assay/mz/rt windows (per-dataset; see body)
#   3. extract_eics(...) → EICs
#   4. process_eics(...) → narrow chromatograms + metrics (FWHM flag)
#   5. build_sample_table(...) → per-file table with detection_tier + log2_ratio
#   6. run_ms1_for_dataset(...) wraps 1, 3, 4, 5 for the cross-dataset loop

# ==========================================================================
# ---- 0. peak-picker core -------------------------------------------------
# Raw (rt, it) of one file + an acceptance RT window -> one integrated peak (or
# NULL), plus the accept gate and the MS1-apex window builder. Pure numerics on
# vectors (no Chromatograms/backend dependency); the pipeline stages below feed it
# each EIC's (rtime, intensity). Design + every threshold: 03-all-datasets.qmd and
# the peak-picker-redesign memory. In brief: new_win() snaps each confirmed MS2 hit
# to its MS1 apex and brackets the robust consensus (+/- APEX_PAD); pick() integrates
# each in-window local max valley-to-valley (merging <PROX_S-from-apex pairs + shallow
# saddles, stopping at deep valleys / dead gaps), caps the tail at +/- MAXW_FWHM*FWHM,
# and selects the dominant peak by AUC weighted (RP) toward the expected RT so an
# on-RT peak beats a larger-area interferent (the isobaric omeprazole sulfone);
# detected_new() is the chromatography-conditional accept gate.

WIDE_EXTRA <- 50      # s of extra EIC each side of the window so the valley walk has data
HWS        <- 1L      # local-maxima half-window (beat immediate neighbours; suits sparse peaks)
MIN_PTS       <- 3L   # RP: sharp peaks legitimately have few acquisitions
MIN_PTS_HILIC <- 4L   # HILIC: coarse sampling (~2 s/scan) gives a real 5 s peak only 4 points
GAP_FACTOR <- 4       # end the walk if the next scan is > this many median spacings away (baseline)
GAP_MIN_S  <- 15      # ... but never below this absolute width (broad sparse HILIC peaks gap ~8 s)
FWHM_BAND_HILIC <- c(1.3, Inf)  # HILIC: reject sub-1.3 s spikes; point-count + slope gates do the rest
FWHM_BAND_RP    <- c(1, 30)     # RP peaks top out ~24 s; 30 admits them, cuts only 40 s+ matrix junk
MERGE_FRAC <- 0.65    # merge across a saddle only if the valley stays above this frac of the lower top
PROX_S     <- 4       # ... but always merge a maximum < this many s FROM THE APEX (one unresolved pair)
NEAR_S     <- 15      # subordinate-peak filter: look for a taller feature within this many s
NBR_MAX    <- 1.3     # ... reject if a local max that close (outside the peak) is > this x the apex
PROM_MIN   <- 0.25    # return-to-baseline: apex must rise >= this frac above its higher boundary
SNR_MIN_HILIC <- 3    # apex height above its local baseline line, in noise units (HILIC floor)
SNR_MIN_RP    <- 12   # RP floor: every real RP peak clears SNR 20; 12 sits below them, above flat ramps
MIN_INT       <- 500  # absolute apex-intensity floor (both modes): below this a peak is noise
MAXW_FWHM  <- 3       # cap each half-boundary at this many FWHM from the apex (drop a sharp peak's tail)
APEX_PAD   <- 5       # pad to LOCATE the apex window (new_win); also the tail-rejection tolerance
RT_SIGMA   <- 4       # RP selection: gaussian down-weight of AUC by distance from the expected RT
SNAP_W     <- 20      # new_win: snap a confirmed MS2 rt to the most intense MS1 apex within +/- this
CLUST_TOL  <- 12      # new_win: keep snapped apexes within this of their median (robust consensus)

# Walk one flank from the apex to the boundary. Descend to a valley; keep going past it (merge the
# shoulder max beyond) when that shoulder sits < PROX_S from the APEX (a close pair, whatever the
# valley depth) OR the saddle is shallow (valley > MERGE_FRAC of the shoulder). Stop at the baseline
# floor or a dead gap (the EIC dropped out -> the peak ended; do not glue on a far point). Proximity
# is from the apex, not the last-merged max, so it cannot chain across a dense trace to a distant peak.
walk_bound <- function(it, ai, dir, floor_lv, apex, rt = NULL, gap = Inf, clust = 0) {
  n <- length(it); i <- ai
  jumped <- function(a, b) !is.null(rt) && abs(rt[b] - rt[a]) > gap
  repeat {
    j <- i + dir
    if (j < 1L || j > n) return(i)
    if (jumped(i, j)) return(i)
    if (it[j] <= it[i]) {
      i <- j; if (it[i] <= floor_lv) return(i)
    } else {
      if (it[i] <= floor_lv) return(i)
      k <- j; while (k + dir >= 1L && k + dir <= n && it[k + dir] >= it[k] && !jumped(k, k + dir)) k <- k + dir
      near <- clust > 0 && !is.null(rt) && abs(rt[k] - rt[ai]) < clust
      if (it[k] > floor_lv && (near || it[i] > MERGE_FRAC * min(it[k], apex))) i <- k
      else return(i)
    }
  }
}

# One file's raw (rt, it) + acceptance window -> one peak (dominant by RT-weighted AUC) or NULL.
pick <- function(rt, it, acc_lo, acc_hi, hilic = FALSE) {
  o0 <- order(rt); rt <- rt[o0]; it <- it[o0]
  gap_thr <- max(GAP_FACTOR * stats::median(diff(rt)), GAP_MIN_S)
  # baseline & noise from the FULL window (NA/below-detection scans read as 0) so the estimate
  # includes the real baseline; detected-points-only inflates noise on a peak-dominated sparse EIC.
  it0 <- it; it0[!is.finite(it0) | it0 < 0] <- 0
  baseline <- as.numeric(stats::quantile(it0, 0.10))
  bpts  <- it0[it0 <= stats::median(it0)]
  noise0 <- if (length(bpts) >= 3) stats::mad(bpts, constant = 1.4826) else NA_real_
  noise_at <- function(a) if (!is.finite(noise0) || noise0 <= 0) 0.01 * max(a - baseline, 1) else noise0
  ok <- is.finite(it) & it > 0
  rt <- rt[ok]; it <- it[ok]
  if (length(rt) < 3L) return(NULL)
  lm <- MsCoreUtils::localMaxima(it, hws = HWS)
  cand <- which(lm & rt >= acc_lo & rt <= acc_hi)
  if (!length(cand)) return(NULL)
  build_peak <- function(ai, resolved = FALSE) {
    # resolve the candidate to the true apex of its region (a local max can sit on a bigger tail)
    if (!resolved) {
      fl0 <- baseline + 3 * noise_at(it[ai])
      i0  <- walk_bound(it, ai, -1L, fl0, it[ai], rt, gap_thr, PROX_S) : walk_bound(it, ai, +1L, fl0, it[ai], rt, gap_thr, PROX_S)
      tai <- i0[which.max(it[i0])]
      if (tai != ai) return(build_peak(tai, TRUE))
    }
    apex <- it[ai]; a_rt <- rt[ai]
    noise <- noise_at(apex)
    floor_lv <- baseline + 3 * noise
    li <- walk_bound(it, ai, -1L, floor_lv, apex, rt, gap_thr, PROX_S)
    ri <- walk_bound(it, ai, +1L, floor_lv, apex, rt, gap_thr, PROX_S)
    # full-peak FWHM: OUTERMOST half-max crossings within [li, ri] (a merged cluster keeps its width)
    half <- (apex + baseline) / 2
    seg <- li:ri; ab <- seg[it[seg] > half]
    fwhm <- if (!length(ab)) NA_real_ else {
      lo <- min(ab); hi <- max(ab)
      xl <- if (lo > li) approx(it[c(lo - 1, lo)], rt[c(lo - 1, lo)], half)$y else NA_real_
      xr <- if (hi < ri) approx(it[c(hi, hi + 1)], rt[c(hi, hi + 1)], half)$y else NA_real_
      if (!is.na(xl) && !is.na(xr)) xr - xl
      else if (!is.na(xl)) 2 * (a_rt - xl)
      else if (!is.na(xr)) 2 * (xr - a_rt)
      else NA_real_
    }
    if (is.finite(fwhm) && fwhm > 0) {
      while (li < ai && rt[li] < a_rt - MAXW_FWHM * fwhm) li <- li + 1
      while (ri > ai && rt[ri] > a_rt + MAXW_FWHM * fwhm) ri <- ri - 1
    }
    idx <- li:ri; npts <- length(idx)
    if (npts < 3L) return(NULL)
    # reject a tail: integrated region's max outside the window -> flank of a bigger peak elsewhere
    tmax_rt <- rt[idx[which.max(it[idx])]]
    if (tmax_rt < acc_lo - APEX_PAD || tmax_rt > acc_hi + APEX_PAD) return(NULL)
    bl  <- if (rt[ri] > rt[li]) it[li] + (it[ri] - it[li]) * (rt[idx] - rt[li]) / (rt[ri] - rt[li])
           else rep(min(it[li], it[ri]), npts)
    ic  <- pmax(it[idx] - bl, 0)
    auc <- sum(diff(rt[idx]) * (head(ic, -1) + tail(ic, -1)) / 2)
    apex_base <- bl[ai - li + 1L]
    nb <- lm & abs(rt - a_rt) <= NEAR_S & !(rt >= rt[li] & rt <= rt[ri])
    nbr_ratio <- if (any(nb)) max(it[nb]) / apex else 0
    list(apex_rt = a_rt, apex = apex, lb = rt[li], rb = rt[ri], idx = idx,
         lb_int = it[li], rb_int = it[ri], npts = npts, auc = auc, fwhm = fwhm,
         prominence = apex - max(it[li], it[ri]), pkwidth = rt[ri] - rt[li], nbr_ratio = nbr_ratio,
         snr = (apex - apex_base) / noise, rt = rt, it = it)
  }
  peaks <- Filter(Negate(is.null), lapply(cand, build_peak))
  if (!length(peaks)) return(NULL)
  # RP: weight AUC by a gaussian in distance from the expected RT (window centre) so the on-RT peak
  # beats a larger-area interferent elsewhere. HILIC: pure area-select (its windows were audited).
  auc <- vapply(peaks, function(p) p$auc, numeric(1))
  score <- if (hilic) auc else
    auc * exp(-0.5 * ((vapply(peaks, function(p) p$apex_rt, numeric(1)) - (acc_lo + acc_hi) / 2) / RT_SIGMA)^2)
  peaks[[which.max(score)]]
}

# Chromatography-conditional accept gate. Shared: point count, intensity floor, SNR, FWHM band,
# slope (FWHM within base width). HILIC also: return-to-baseline (prominence) + no much-taller
# neighbour (both audited on HILIC). RP drops those two (they killed sharp truncated-flank peaks
# and peaks beside the sulfone) but still demands a 3-point peak clear its own baseline.
detected_new <- function(p, hilic) {
  if (is.null(p) || !is.finite(p$auc)) return(FALSE)
  min_pts <- if (hilic) MIN_PTS_HILIC else MIN_PTS
  fw      <- if (hilic) FWHM_BAND_HILIC else FWHM_BAND_RP
  snr_min <- if (hilic) SNR_MIN_HILIC  else SNR_MIN_RP
  ok <- p$npts >= min_pts && p$apex >= MIN_INT &&
    is.finite(p$snr) && p$snr >= snr_min &&
    is.finite(p$fwhm) && p$fwhm >= fw[1] && p$fwhm <= fw[2] &&
    p$fwhm <= 1.8 * p$pkwidth
  if (!ok) return(FALSE)
  baseline_ok <- p$prominence >= PROM_MIN * p$apex
  if (hilic) baseline_ok && p$nbr_ratio <= NBR_MAX
  else       p$npts >= 4L || baseline_ok
}

# snap one confirmed MS2 rt to the most intense MS1 apex within +/- SNAP_W of it
snap_one <- function(rt, it, center) {
  ok <- is.finite(rt) & is.finite(it) & it > 0; rt <- rt[ok]; it <- it[ok]
  if (length(rt) < 3) return(NA_real_)
  o <- order(rt); rt <- rt[o]; it <- it[o]
  lm <- MsCoreUtils::localMaxima(it, hws = 1L)
  cand <- which(lm & abs(rt - center) <= SNAP_W); if (!length(cand)) return(NA_real_)
  rt[cand[which.max(it[cand])]]
}

# MS1-apex acceptance window: consensus of the snapped confirmed apexes +/- APEX_PAD. `cache_role`
# is list(files = basenames, traces = list(list(rt=, it=))); `bare` normalises a filename to the
# match key; `old_win` is the fallback. `conf` is the confirmed-hit table (== conf_final, cols
# dataset_id/role/rt/source_file); `ov` the manual apex overrides (rt_apex_overrides.csv).
new_win <- function(ds, role_full, cache_role, bare, old_win, conf, ov) {
  hh <- conf[conf$dataset_id == ds & conf$role == role_full, ]
  if (!nrow(hh)) return(old_win)
  fk <- bare(cache_role$files); ms2 <- suppressWarnings(as.numeric(hh$rt))
  ovr <- ov[ov$dataset_id == ds & ov$role == role_full, ]
  apex <- vapply(seq_len(nrow(hh)), function(i) {
    j <- match(bare(hh$source_file[i]), fk); if (is.na(j)) return(NA_real_)
    k <- which(abs(ms2[i] - ovr$ms2_rt) <= 2)
    if (length(k)) return(ovr$apex_rt[k[1]])
    snap_one(cache_role$traces[[j]]$rt, cache_role$traces[[j]]$it, ms2[i])
  }, numeric(1))
  a <- apex[is.finite(apex)]; if (!length(a)) return(old_win)
  cons <- stats::median(a); inl <- a[abs(a - cons) <= CLUST_TOL]
  if (!length(inl)) return(old_win)
  pad <- if (length(inl) >= 2) APEX_PAD else 2 * APEX_PAD
  c(min(inl) - pad, max(inl) + pad)
}

# Manual reviewer drop-list, keyed by (dataset, role, file) on norm_file_key -- a stable key that
# survives re-extraction and file re-ordering. `drop_all` = ms1_drop_list.csv (cols dataset/role/file).
is_dropped <- function(dataset, role_full, file, drop_all) {
  if (!nrow(drop_all)) return(FALSE)
  dl <- drop_all[drop_all$dataset == dataset & drop_all$role == role_full, ]
  nrow(dl) > 0 && norm_file_key(file) %in% norm_file_key(dl$file)
}

# ---- 1. spectrum loading (repo dispatch) ---------------------------------

load_spectra_ms1 <- function(deposit_id, files, assay = character()) {
  wanted <- basename(files)
  pat    <- paste0("(", paste(gsub("\\.", "\\\\.", wanted),
                              collapse = "|"), ")$")
  be <- if (startsWith(deposit_id, "MTBLS")) {
          # Pass assayName so the correct assay's files load. MsBackendMetaboLights
          # (>= 1.7.4, gabri) bakes the assay index into the cache filename, so a
          # deposit's pos/neg files with identical basenames no longer collide.
          MsBackendMetaboLights::mtbls_sync_data_files(
            mtblsId = deposit_id, assayName = assay, fileName = wanted)
          backendInitialize(MsBackendMetaboLights(),
                            mtblsId   = deposit_id,
                            assayName = assay,
                            fileName  = wanted,
                            offline   = TRUE)
        } else if (startsWith(deposit_id, "ST")) {
          MsBackendMetabolomicsWorkbench::mwb_sync_data_files(
            mwbId = deposit_id, fileName = wanted)
          backendInitialize(MsBackendMetabolomicsWorkbench(),
                            mwbId       = deposit_id,
                            filePattern = pat,
                            offline     = TRUE)
        } else if (startsWith(deposit_id, "MSV")) {
          MsBackendMassIVE::massive_sync_data_files(
            massiveId = deposit_id, fileName = wanted)
          backendInitialize(MsBackendMassIVE(),
                            massiveId   = deposit_id,
                            filePattern = pat,
                            offline     = TRUE)
        } else
          stop("Unknown deposit_id prefix: ", deposit_id)
  filterMsLevel(Spectra(be), 1L)
}

# ---- 2. derive mz/rt windows from confirmed features ---------------------

dominant_adduct <- function(df) {
  tab <- table(df$adduct_inferred[!is.na(df$adduct_inferred)])
  if (!length(tab)) return(NA_character_)
  names(tab)[which.max(tab)]
}

# Pick the assay with the most confirmed (parent+metabolite) features, then
# derive m/z and RT windows from the dominant adduct per role. Returns a list
# ready for extract_eics().
# rt_pad_s = 5 is not a comfort default: it is load-bearing. 5-hydroxyomeprazole
# (CYP2C19) and omeprazole sulfone (CYP3A4) are exact isomers (C17H19N3O4S,
# m/z 362.1169), so MS1 cannot tell them apart and only retention time separates
# them. In MSV000084008 the confirmed 5-OH span is 3.4 s wide (154.9-158.3) and
# the sulfone sits at ~181 s. A wide pad lets the sulfone into the window, and
# peakBoundary() then integrates it instead (it is the larger peak in a poor
# CYP2C19 metaboliser). Validated against the source study's published
# PharmKinetics_CYP2C19: pad 3-20 s gives r = +0.55 to +0.62 (p < 0.005), pad 30 s
# gives r = -0.01 (p = 0.96), i.e. the signal is destroyed exactly when the pad
# reaches the sulfone. Do not widen this without re-running that check.
derive_extraction_config <- function(dataset_id, conf_final,
                                     rt_pad_s = 5, bio_files = NULL) {
  pc <- conf_final[conf_final$dataset_id == dataset_id, ]
  if (!nrow(pc)) stop("No confirmed features for ", dataset_id)

  # Membership rule: the retention-time anchor must come only from confirmed
  # hits on files we actually extract. A MASST/SIRIUS hit can live on a file that
  # is not in this deposit's extraction set (e.g. MSV000082493's confirmations
  # sit on blood files that belong to the re-deposit MSV000084008, while 082493
  # itself contributes feces/skin/urine). Anchoring on those imports a foreign
  # matrix's retention time. Filter here so the window is always derived from the
  # same population it is applied to. A dataset left with no parent or no
  # metabolite anchor after filtering has no valid in-set confirmation and must
  # be re-confirmed (see the biofluid-stratified cap in the confirmation qmds),
  # so we stop loudly rather than fall back to foreign files.
  if (!is.null(bio_files)) {
    # Match confirmed-hit files to the extraction set on a normalised key that
    # drops our accession prefix (MsBackend caches MassIVE files as MSV000..._<name>;
    # we prefix MetaboLights/Workbench downloads with the accession), the
    # extension (.mzML vs .mzXML re-deposits), and case. Same rule as
    # norm_file_key() in R/metadata.R.
    stem <- function(x) tolower(sub("[.][^.]*$", "",
              sub("^(MSV[0-9]+|MTBLS[0-9]+|ST[0-9]+)[_-]", "",
                  basename(as.character(x)))))
    in_set <- stem(pc$source_file) %in% stem(bio_files)
    n_par <- sum(in_set & pc$role == "parent")
    n_met <- sum(in_set & pc$role == "metabolite")
    if (n_par == 0 || n_met == 0)
      stop("Membership filter leaves ", dataset_id, " with par=", n_par,
           " met=", n_met, " in-set confirmations; needs re-confirmation on ",
           "its own extraction set before extraction.")
    pc <- pc[in_set, ]
  }

  assay_tab <- pc |>
    dplyr::group_by(assay) |>
    dplyr::summarise(n_par   = sum(role == "parent"),
                     n_met   = sum(role == "metabolite"),
                     .groups = "drop") |>
    dplyr::filter(n_par > 0, n_met > 0) |>
    dplyr::mutate(n_total = n_par + n_met) |>
    dplyr::arrange(desc(n_total))
  if (!nrow(assay_tab))
    stop("No assay with both parent and metabolite confirmed in ", dataset_id)

  assay <- assay_tab$assay[1]
  ac    <- if (is.na(assay)) pc[is.na(pc$assay), ]
           else              pc[!is.na(pc$assay) & pc$assay == assay, ]
  ad_p  <- dominant_adduct(ac[ac$role == "parent",     ])
  ad_m  <- dominant_adduct(ac[ac$role == "metabolite", ])

  mz_p <- median(ac$mz[ac$role == "parent"     & ac$adduct_inferred == ad_p],
                 na.rm = TRUE)
  mz_m <- median(ac$mz[ac$role == "metabolite" & ac$adduct_inferred == ad_m],
                 na.rm = TRUE)
  # Per-dataset RT window: span this dataset's confirmed RTs (min..max) and pad
  # each side by rt_pad_s. The confirmed span already absorbs cross-file drift
  # (it is the spread of MS2 detections across many files), so the pad only has
  # to cover drift in files that carry no confirmed hit. Keep it small: see the
  # isomer note above for what a generous pad costs.
  rt_p <- ac$rt[ac$role == "parent"     & ac$adduct_inferred == ad_p]
  rt_m <- ac$rt[ac$role == "metabolite" & ac$adduct_inferred == ad_m]
  rt_p <- rt_p[is.finite(rt_p)]; rt_m <- rt_m[is.finite(rt_m)]

  list(
    assay      = assay,
    adduct_par = ad_p, adduct_met = ad_m,
    mz_par     = mz_p, mz_met     = mz_m,
    rt_par_win = c(min(rt_p) - rt_pad_s, max(rt_p) + rt_pad_s),
    rt_met_win = c(min(rt_m) - rt_pad_s, max(rt_m) + rt_pad_s),
    confirmed_files_parent     = unique(basename(
      ac$source_file[ac$role == "parent"     & !is.na(ac$source_file)])),
    confirmed_files_metabolite = unique(basename(
      ac$source_file[ac$role == "metabolite" & !is.na(ac$source_file)]))
  )
}

# ---- 3. EIC extraction ---------------------------------------------------

# Padding kept in chr_all beyond the detection windows. Must stay >= the widest
# re-extraction any consumer does against chr_all: dataset-report.qmd inspects
# "neither" files at +/- 60 s (WIDE_PAD_S), so 75 s leaves headroom. Narrowing this
# further shrinks chr_all more but would silently return empty EICs there.
CHR_ALL_PAD_S <- 75

extract_eics <- function(ms1, mz_par, mz_met,
                         rt_par_win, rt_met_win, ppm = 20) {
  # Restrict the spectra to the RT/mz neighbourhood before building chr_all.
  # Verified on MSV000084008 (40 files, 29 detections): identical files and a max
  # relative AUC difference of 0, with chr_all much smaller and extract_eics ~16%
  # faster. The dominant cost is now the per-file pick() loop, not the extraction.
  rt_lo <- min(rt_par_win[1], rt_met_win[1]) - CHR_ALL_PAD_S
  rt_hi <- max(rt_par_win[2], rt_met_win[2]) + CHR_ALL_PAD_S
  mz_lo <- min(mz_par, mz_met) * (1 - 50e-6)
  mz_hi <- max(mz_par, mz_met) * (1 + 50e-6)
  ms1   <- filterMzRange(filterRt(ms1, c(rt_lo, rt_hi)), c(mz_lo, mz_hi))

  all_orig <- unique(dataOrigin(ms1))
  peak_tbl <- rbind(
    # Extract WIDE_EXTRA s past each acceptance window so pick()'s valley-to-valley walk always has
    # data (broad peaks run far past the apex); new_win/pick then only accept apexes in the window.
    data.frame(role = "parent",     msLevel = 1L, dataOrigin = all_orig,
               rtMin = rt_par_win[1] - WIDE_EXTRA, rtMax = rt_par_win[2] + WIDE_EXTRA,
               mzMin = mz_par * (1 - ppm * 1e-6),
               mzMax = mz_par * (1 + ppm * 1e-6),
               stringsAsFactors = FALSE),
    data.frame(role = "metabolite", msLevel = 1L, dataOrigin = all_orig,
               rtMin = rt_met_win[1] - WIDE_EXTRA, rtMax = rt_met_win[2] + WIDE_EXTRA,
               mzMin = mz_met  * (1 - ppm * 1e-6),
               mzMax = mz_met  * (1 + ppm * 1e-6),
               stringsAsFactors = FALSE)
  )
  chr_all  <- Chromatograms(ms1)
  # chromExtract yields a lazily Spectra-backed object, so every downstream
  # peaksData()/intensity() re-materialises it from the spectra — count_peaks() and the
  # pick() trace-build each paid that (~15% of extract compute on MSV000084008). Pull the
  # targeted EICs into memory once with setBackend(): subsequent access (count_peaks, pick,
  # dataset-report plots) then reads from RAM. chr_all is left Spectra-backed — it is only
  # re-chromExtract'd (reports/panels) and stays small in the cache.
  chr_eics <- setBackend(
    chromExtract(chr_all, peak_tbl, by = c("msLevel", "dataOrigin")),
    ChromBackendMemory())
  list(
    chr_all = chr_all,
    chr_par = chr_eics[chromData(chr_eics)$role == "parent"],
    chr_met = chr_eics[chromData(chr_eics)$role == "metabolite"]
  )
}

# ---- 4. peak pick (per file) + metrics + FWHM flag -----------------------

# For each file: build the MS1-apex acceptance window (new_win), pick the dominant peak (pick),
# gate it (detected_new) and apply the manual drop-list. Emits metrics_* with the SAME schema as
# before (auc/fwhm/prominence/pkwidth/snr) PLUS a `detected` flag, the apex RT and point count, so
# build_sample_table and 05 read one authoritative detection instead of re-deriving a shape rule.
# count_peaks and the (wide) per-file EICs are kept as chr_*_narrow for dataset-report diagnostics.
ROLE_COLS <- c("auc", "fwhm", "prominence", "pkwidth", "snr", "apex_rt", "npts", "detected")

process_eics <- function(chr_par, chr_met, chr_all, dataset_id, cfg, conf, ov, drop_all,
                         is_hilic, files_base) {
  peaks_par <- count_peaks(chr_par)
  peaks_met <- count_peaks(chr_met)
  bare <- function(x) norm_file_key(to_bare_name(basename(as.character(x)), files_base))

  role_metrics <- function(chr, role, role_full, old_win) {
    files  <- basename(dataOrigin(chr))
    traces <- lapply(Chromatograms::peaksData(chr),
                     function(pd) list(rt = pd[, "rtime"], it = pd[, "intensity"]))
    win    <- new_win(dataset_id, role_full, list(files = files, traces = traces),
                      bare, old_win, conf, ov)
    picks  <- lapply(traces, function(tr) pick(tr$rt, tr$it, win[1], win[2], is_hilic))
    keep   <- which(!vapply(picks, is.null, logical(1)))
    cols   <- paste0(role, "_", ROLE_COLS)
    if (!length(keep)) {
      m <- data.frame(file = character(0), stringsAsFactors = FALSE)
      for (cc in cols) m[[cc]] <- if (endsWith(cc, "detected")) logical(0) else numeric(0)
      return(m)
    }
    pk <- picks[keep]; fl <- files[keep]
    det <- vapply(seq_along(pk), function(i)
      detected_new(pk[[i]], is_hilic) &&
      !is_dropped(dataset_id, role_full, fl[i], drop_all), logical(1))
    m <- data.frame(
      file       = fl,
      auc        = vapply(pk, function(p) p$auc,        numeric(1)),
      fwhm       = vapply(pk, function(p) p$fwhm,       numeric(1)),
      prominence = vapply(pk, function(p) p$prominence, numeric(1)),
      pkwidth    = vapply(pk, function(p) p$pkwidth,    numeric(1)),
      snr        = vapply(pk, function(p) p$snr,        numeric(1)),
      apex_rt    = vapply(pk, function(p) p$apex_rt,    numeric(1)),
      npts       = vapply(pk, function(p) as.integer(p$npts), integer(1)),
      detected   = det,
      stringsAsFactors = FALSE)
    stats::setNames(m, c("file", cols))
  }

  metrics_par <- role_metrics(chr_par, "par", "parent",     cfg$rt_par_win)
  metrics_met <- role_metrics(chr_met, "met", "metabolite", cfg$rt_met_win)
  metrics_par$par_fwhm_outlier <- flag_fwhm(metrics_par, "par_fwhm")
  metrics_met$met_fwhm_outlier <- flag_fwhm(metrics_met, "met_fwhm")

  list(peaks_par = peaks_par, peaks_met = peaks_met,
       chr_par_narrow = chr_par, chr_met_narrow = chr_met,  # wide EICs serve dataset-report plots
       metrics_par = metrics_par, metrics_met = metrics_met)
}

# ---- 5. sample-table assembly --------------------------------------------

# NOTE: to_bare_name() now lives in R/metadata.R, next to norm_file_key(), which
# calls it. Keeping it here would make metadata.R depend on ms1.R being sourced
# first, and several scripts source metadata.R alone.

build_sample_table <- function(metrics_par, metrics_met, files_base, meta,
                               dedup = TRUE) {
  # dedup = TRUE  -> one row per subject (default; 05 clinical-group analyses,
  #                  avoids pseudoreplication).
  # dedup = FALSE -> one row per file, subject_id retained, so genotype /
  #                  longitudinal / specimen analyses dedup themselves as needed
  #                  (see 05's per-file table and 05-cyp-validation.qmd).
  # Reconcile prefixed cache filenames to the bare bio_files names so the joins
  # below match. Without this every metrics row misses and the dataset collapses
  # to "neither".
  metrics_par$file <- to_bare_name(metrics_par$file, files_base)
  metrics_met$file <- to_bare_name(metrics_met$file, files_base)
  sample_tbl <- merge(metrics_par, metrics_met, by = "file", all = TRUE)
  sample_tbl <- merge(data.frame(file = files_base, stringsAsFactors = FALSE),
                      sample_tbl, by = "file", all.x = TRUE)

  # Join the metadata on the NORMALISED key, not the raw name. A metadata source
  # may spell the extension differently from the deposit: MSV000084008's submitter
  # TSV lists "DM001_9_000_RD11_01_43875.mzXML" while the deposit holds the same
  # file as ".mzML". Merging on the literal string matched 0 of 461 rows there, so
  # subject_id, body_site and timepoint all silently became NA for the entire
  # plasma cohort. norm_file_key drops the extension (and the accession), and
  # matches 439.
  sample_tbl$.k <- norm_file_key(sample_tbl$file)
  meta          <- meta[!duplicated(norm_file_key(meta$file)), , drop = FALSE]
  meta$.k       <- norm_file_key(meta$file)
  meta$file     <- NULL                     # keep the deposit's spelling
  sample_tbl    <- merge(sample_tbl, meta, by = ".k", all.x = TRUE)
  sample_tbl$.k <- NULL

  # Detection is decided once, in the picker (detected_new + the manual drop-list); metrics_*
  # carry it as par_detected / met_detected. Files with no picked peak merge in as NA -> FALSE.
  # The gate itself (chromatography-conditional) is documented in 03-all-datasets.qmd.
  par_detected <- sample_tbl$par_detected %in% TRUE
  met_detected <- sample_tbl$met_detected %in% TRUE

  sample_tbl$detection_tier <-
    ifelse( par_detected &  met_detected, "both",
    ifelse( par_detected & !met_detected, "parent_only",
    ifelse(!par_detected &  met_detected, "metab_only",
                                          "neither")))

  sample_tbl$log2_ratio <- ifelse(
    sample_tbl$detection_tier == "both",
    log2(sample_tbl$par_auc / sample_tbl$met_auc), NA_real_)

  # Per-file table (subject_id retained) when dedup is off.
  if (!dedup) return(as.data.frame(sample_tbl))

  # Dedup to one file per subject by highest combined AUC. Fall back to
  # sample_name then file when subject_id is NA, else all NA rows collapse to
  # one group.
  sample_tbl$.dedup_key <- dplyr::coalesce(
    as.character(sample_tbl$subject_id),
    as.character(sample_tbl$sample_name),
    as.character(sample_tbl$file))

  sample_tbl |>
    dplyr::group_by(.dedup_key) |>
    dplyr::slice_max(
      order_by = dplyr::coalesce(par_auc, 0) + dplyr::coalesce(met_auc, 0),
      n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-.dedup_key) |>
    as.data.frame()
}

# ---- 6. cross-dataset wrapper --------------------------------------------

# Cross-dataset entry point. Resolves deposit_id, assay, m/z and RT windows
# from the artifacts; the caller supplies only the dataset_id and three artifacts.
run_ms1_for_dataset <- function(dataset_id,
                                conf_final,
                                bio_files_per_ds,
                                ds_meta_tbl,
                                ppm           = 20,
                                rt_pad_s      = 5,
                                verbose       = TRUE) {
  row <- ds_meta_tbl[ds_meta_tbl$dataset_id == dataset_id, ]
  if (!nrow(row))
    stop("No row for ", dataset_id, " in dataset_metadata_table")
  deposit_id <- row$deposit_id
  assay      <- row$assay

  files <- bio_files_per_ds[[dataset_id]]
  cfg   <- derive_extraction_config(dataset_id, conf_final,
                                    rt_pad_s = rt_pad_s, bio_files = files)
  meta  <- load_metadata(dataset_id)
  # Picker inputs: HILIC vs RP gates, manual apex overrides, and the reviewer drop-list.
  is_hilic <- { chrom <- row$chromatography; !is.na(chrom) && toupper(chrom) == "HILIC" }
  ov <- if (file.exists("curation/rt_apex_overrides.csv"))
    read.csv("curation/rt_apex_overrides.csv", stringsAsFactors = FALSE) else
    data.frame(dataset_id = character(), role = character(), ms2_rt = numeric(), apex_rt = numeric())
  drop_all <- if (file.exists("curation/ms1_drop_list.csv"))
    read.csv("curation/ms1_drop_list.csv", stringsAsFactors = FALSE) else
    data.frame(dataset = character(), role = character(), file = character(), reason = character())

  if (verbose) message("[", dataset_id, "] loading ", length(files),
                       " MS1 files …")
  # For MetaboLights, pass the assay so same-basename pos/neg modes in separate
  # assays are disambiguated at fetch time (ignored for non-MTBLS deposits).
  ms1 <- load_spectra_ms1(deposit_id, files, assay = assay)

  if (verbose) message("[", dataset_id, "] extracting EICs …")
  eics <- extract_eics(ms1, cfg$mz_par, cfg$mz_met,
                       cfg$rt_par_win, cfg$rt_met_win, ppm)

  if (verbose) message("[", dataset_id, "] picking peaks + metrics …")
  proc <- process_eics(eics$chr_par, eics$chr_met, eics$chr_all, dataset_id, cfg,
                       conf_final, ov, drop_all, is_hilic, basename(files))

  par_files_all <- basename(dataOrigin(eics$chr_par))
  met_files_all <- basename(dataOrigin(eics$chr_met))
  if (length(cfg$confirmed_files_parent))
    check_confirmed(par_files_all, proc$peaks_par,
                    cfg$confirmed_files_parent, "Parent    ")
  if (length(cfg$confirmed_files_metabolite))
    check_confirmed(met_files_all, proc$peaks_met,
                    cfg$confirmed_files_metabolite, "Metabolite")

  sample_tbl <- build_sample_table(proc$metrics_par, proc$metrics_met,
                                    basename(files), meta)

  list(
    dataset_id     = dataset_id,
    deposit_id     = deposit_id,
    assay          = assay,
    cfg            = cfg,
    chr_all        = eics$chr_all,
    chr_par        = eics$chr_par,
    chr_met        = eics$chr_met,
    chr_par_narrow = proc$chr_par_narrow,
    chr_met_narrow = proc$chr_met_narrow,
    peaks_par      = proc$peaks_par,
    peaks_met      = proc$peaks_met,
    metrics_par    = proc$metrics_par,
    metrics_met    = proc$metrics_met,
    sample_tbl     = sample_tbl
  )
}

# ===== EIC extraction (was eic.R) =====

# R/eic.R — EIC extraction helpers

# Compute a shared y-axis limit across all EICs in a Chromatograms object.
# The 5% headroom prevents the apex from touching the plot border.
eic_ylim <- function(chr) {
  ints <- unlist(intensity(chr), use.names = FALSE)
  c(0, max(ints, na.rm = TRUE) * 1.05)
}

# ===== peak detection + metrics (was peaks.R) =====

# R/peaks.R — peak detection and EIC quality metrics

# Count chromatographic peaks in each EIC using a local-maxima detector.
# hws: half-window size in scans (a peak must be the maximum within ±hws
# neighbours). pct_max: minimum height as a fraction of the trace maximum
# (filters sub-noise wiggles).
count_peaks <- function(chr, hws = 3L, pct_max = 0.1) {
  ints <- intensity(chr)
  vapply(ints, function(x) {
    mx <- max(x, na.rm = TRUE)
    if (!is.finite(mx)) return(0L)
    x_s <- x[!is.na(x)]
    lm  <- MsCoreUtils::localMaxima(x_s, hws = hws)
    sum(lm & x_s >= pct_max * mx)
  }, integer(1L))
}

# Flag files whose FWHM deviates more than 3 × MAD from the dataset median.
# NA inputs (files with no detected peak) are treated as non-outliers so they
# do not propagate into downstream sum() calls.
flag_fwhm <- function(df, col) {
  v   <- df[[col]]
  med <- median(v, na.rm = TRUE)
  mad <- stats::mad(v, na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) return(rep(FALSE, length(v)))
  flagged <- abs(v - med) > 3 * mad
  flagged[is.na(flagged)] <- FALSE
  flagged
}

# Report whether SIRIUS-confirmed files produced a detected peak. Confirmed
# files are positive controls: a miss is a settings problem (ppm, RT window,
# adduct), not biology. Extension-agnostic match — MassIVE often hosts the same
# scan data as both .mzML and .mzXML, and submitter vs Pan-ReDU TSV can disagree.
check_confirmed <- function(files, n_peaks, conf_files, role) {
  # Normalise both sides the way the rest of the pipeline does (drop our accession
  # prefix, extension, case) so a naming difference is never mistaken for a miss.
  key       <- function(x) tolower(sub("[.][^.]*$", "",
                 sub("^(MSV[0-9]+|MTBLS[0-9]+|ST[0-9]+)[_-]", "", basename(x))))
  idx       <- which(key(files) %in% key(conf_files))
  cat(role, "- confirmed files loaded:", length(idx),
      "| with >=1 peak:", sum(n_peaks[idx] > 0), "\n")

  missed <- files[idx[n_peaks[idx] == 0]]
  if (length(missed))
    cat("  WARNING missed:", paste(missed, collapse = ", "), "\n")
}

# ===== within-study stats (was within_stats.R) =====

# R/within_stats.R — within-dataset descriptive + nonparametric tests on a
# sample_tbl from build_sample_table(). Returns a list usable as a per-dataset
# print-out (04) and as a row for the cross-dataset aggregation (05).
#
# Three tests, each guarded for n / group-count:
#   - Fisher exact on parent presence (yes/no) x group
#   - Fisher exact on metab  presence (yes/no) x group
#   - Kruskal-Wallis on log2_ratio x group (both-tier only, FWHM-cleaned)
# Censored cases (parent_only / metab_only) count as "present" in the Fisher
# tests; they can't enter Kruskal-Wallis without a Tobit-style model (future work).

within_dataset_stats <- function(sample_tbl, group_col = "health_status") {
  if (!group_col %in% names(sample_tbl)) {
    warning(group_col, " not in sample_tbl — skipping within-dataset stats")
    return(NULL)
  }
  g <- sample_tbl[[group_col]]
  n_groups <- length(unique(g[!is.na(g)]))

  parent_present <- !is.na(sample_tbl$par_auc)
  metab_present  <- !is.na(sample_tbl$met_auc)

  safe_fisher <- function(present) {
    if (n_groups < 2) return(NULL)
    tab <- table(present = present, group = g)
    if (any(dim(tab) < 2)) return(NULL)
    tryCatch(fisher.test(tab, simulate.p.value = TRUE),
             error = function(e) NULL)
  }

  both <- sample_tbl[sample_tbl$detection_tier == "both" &
                     !sample_tbl$par_fwhm_outlier %in% TRUE &
                     !sample_tbl$met_fwhm_outlier %in% TRUE, ]
  both_g <- both[[group_col]]
  both_n_groups <- length(unique(both_g[!is.na(both_g)]))

  kruskal_ratio <- if (both_n_groups >= 2 && nrow(both) >= 5)
    tryCatch(kruskal.test(both$log2_ratio ~ both_g),
             error = function(e) NULL) else NULL

  per_group_ratio <- if (nrow(both) > 0)
    dplyr::group_by(both, .data[[group_col]]) |>
    dplyr::summarise(
      n          = dplyr::n(),
      median_l2r = median(log2_ratio, na.rm = TRUE),
      iqr_l2r    = stats::IQR(log2_ratio, na.rm = TRUE),
      .groups    = "drop") else NULL

  list(
    group_col       = group_col,
    detection_table = table(detection_tier = sample_tbl$detection_tier,
                            group          = g,
                            useNA          = "ifany"),
    fisher_parent   = safe_fisher(parent_present),
    fisher_metab    = safe_fisher(metab_present),
    kruskal_ratio   = kruskal_ratio,
    per_group_ratio = per_group_ratio
  )
}

# Turn a within_dataset_stats() result into a short plain-English readout.
# Returns a single markdown string; intended for a `results = "asis"` chunk.
interpret_within_stats <- function(ws, alpha = 0.05) {
  if (is.null(ws)) return("_No stratifier available — within-dataset stats skipped._")
  gc <- ws$group_col

  p_par <- if (!is.null(ws$fisher_parent)) ws$fisher_parent$p.value else NA
  p_met <- if (!is.null(ws$fisher_metab))  ws$fisher_metab$p.value  else NA
  p_kw  <- if (!is.null(ws$kruskal_ratio)) ws$kruskal_ratio$p.value else NA

  sig <- function(p) !is.na(p) && p < alpha
  fmt <- function(p) if (is.na(p)) "n/a" else signif(p, 2)

  det_bits <- c()
  if (!is.na(p_par))
    det_bits <- c(det_bits, sprintf("parent detection %s %s (Fisher p = %s)",
      if (sig(p_par)) "**is** associated with" else "is not associated with",
      gc, fmt(p_par)))
  if (!is.na(p_met))
    det_bits <- c(det_bits, sprintf("metabolite detection %s %s (p = %s)",
      if (sig(p_met)) "**is** associated with" else "is not associated with",
      gc, fmt(p_met)))
  det_sentence <- if (length(det_bits))
    paste0("Presence/absence: ", paste(det_bits, collapse = "; "), ".")
  else "Presence/absence: not testable (single group or no variation)."

  ratio_sentence <- if (is.na(p_kw)) {
    "Ratio: not enough `both`-tier samples across ≥2 groups for a Kruskal-Wallis test."
  } else if (sig(p_kw)) {
    sprintf("Ratio: log2(par/met) **differs** by %s among detected samples (Kruskal-Wallis p = %s).",
            gc, fmt(p_kw))
  } else {
    sprintf("Ratio: among detected samples, log2(par/met) does not differ by %s (Kruskal-Wallis p = %s).",
            gc, fmt(p_kw))
  }

  # Combined-pattern note when detection is group-linked but ratio is not
  note <- if (sig(p_par) && !is.na(p_kw) && !sig(p_kw))
    paste0(" Whether the drug is present tracks the clinical grouping, ",
           "but the CYP2C19 ratio in those who carry it does not — ",
           "consistent with prescription pattern (who is on a PPI) driving ",
           "detection, while the metaboliser phenotype is independent of it.")
  else ""

  paste0("**Interpretation.** ", det_sentence, " ", ratio_sentence, note)
}

# Run within_dataset_stats() for every stratifier the sample_tbl supports,
# health_status first (primary clinical group), then the others per
# usable_stratifiers(). Returns a named list keyed by stratifier; empty if none.
within_dataset_stats_multi <- function(sample_tbl) {
  primary <- "health_status"
  cols <- if (primary %in% names(sample_tbl))
            c(primary, usable_stratifiers(sample_tbl, exclude = primary))
          else
            usable_stratifiers(sample_tbl, exclude = character(0))
  if (!length(cols)) return(list())
  out <- vector("list", length(cols))
  names(out) <- cols
  for (gc in cols) out[[gc]] <- within_dataset_stats(sample_tbl, group_col = gc)
  out
}

# Which canonical fields are usable as stratifiers: categorical, present, with
# >= min_levels distinct values each backed by >= min_per_level samples.
# `exclude` drops fields handled elsewhere (health_status by default).
# subject_id / sample_name / file are never stratifiers; age is continuous.
usable_stratifiers <- function(sample_tbl, exclude = "health_status",
                               min_levels = 2L, min_per_level = 3L) {
  cand <- intersect(c("sex", "group", "country", "body_site", "treatment"),
                    names(sample_tbl))
  cand <- setdiff(cand, exclude)
  keep <- character(0)
  for (f in cand) {
    v  <- sample_tbl[[f]]
    v  <- v[!is.na(v) & v != ""]
    tb <- table(v)
    if (length(tb) >= min_levels && sum(tb >= min_per_level) >= min_levels)
      keep <- c(keep, f)
  }
  keep
}

# ===== plots (was plots.R) =====

# R/plots.R — shared plotting functions

# Plot EICs for FWHM-outlier files alongside a reference EIC (the file closest
# to the median FWHM). A genuine outlier looks deformed vs the reference (double
# apex, asymmetric tail, broad shoulder, noisy baseline); if it looks similar,
# the flag is a false alarm and the file should stay in the ratio.
plot_fwhm_outliers <- function(chr, metrics, prefix, role_label) {
  oc        <- paste0(prefix, "_fwhm_outlier")
  fc        <- paste0(prefix, "_fwhm")
  out_files <- metrics$file[metrics[[oc]] %in% TRUE]
  if (!length(out_files) || !length(chr)) {
    cat("No FWHM outliers for ", role_label, " — nothing to plot.\n", sep = "")
    return(invisible())
  }
  files_chr <- basename(dataOrigin(chr))
  med_fwhm  <- median(metrics[[fc]], na.rm = TRUE)
  ref_file  <- metrics$file[which.min(abs(metrics[[fc]] - med_fwhm))]
  to_show   <- unique(c(ref_file, out_files))

  rows <- lapply(seq_along(files_chr), function(i) {
    f <- files_chr[i]
    if (!(f %in% to_show)) return(NULL)
    pd     <- peaksData(chr)[[i]]
    fwhm_v <- metrics[[fc]][match(f, metrics$file)]
    data.frame(
      file      = f,
      rtime     = pd$rtime,
      intensity = pd$intensity,
      panel     = paste0(f, "\nFWHM = ", round(fwhm_v, 1), " s"),
      kind      = ifelse(f == ref_file, "reference (median FWHM)", "outlier"),
      stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || !nrow(df)) return(invisible())

  ggplot(df, aes(x = rtime, y = intensity, colour = kind)) +
    geom_line() +
    facet_wrap(~panel, scales = "free_y") +
    scale_colour_manual(values = c("reference (median FWHM)" = "grey50",
                                   outlier = "#D6604D")) +
    labs(x = "RT (s)", y = "Intensity",
         title  = paste(role_label, "— FWHM outlier EICs"),
         colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position    = "bottom",
          strip.background   = element_rect(fill = "grey92"))
}
