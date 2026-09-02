# R/ms1.R — MS1 EIC → peak → metric pipeline (+ EIC/peaks/stats/plots helpers).
# Stages: load_spectra_ms1 → derive_extraction_config → extract_eics →
# process_eics → build_sample_table; run_ms1_for_dataset wraps them per dataset.

# ==========================================================================
# ---- 0. peak-picker core -------------------------------------------------
# One file's (rt, it) + acceptance window -> one integrated peak (or NULL), the
# accept gate, and the apex-window builder. Pure vector numerics, no backend dep.
# Design + every threshold: 03-all-datasets.qmd and the peak-picker-redesign memory.

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
APEX_PAD   <- 5       # pad to LOCATE the apex window (acceptance_window); also the tail-rejection tolerance
RT_SIGMA   <- 4       # RP selection: gaussian down-weight of AUC by distance from the expected RT
SNAP_W     <- 20      # acceptance_window: snap a confirmed MS2 rt to the most intense MS1 apex within +/- this
CLUST_TOL  <- 12      # acceptance_window: keep snapped apexes within this of their median (robust consensus)
COHERE_S   <- 1       # acceptance_window: extend the window over peaks this close to its edge (cluster cohesion)
CHAN_SEP   <- 6       # below this the two channels are not chromatographically resolved and
                      #   resolve_isobars() cannot tell them apart, so it stands down

# Walk one flank from the apex: descend to a valley, then merge the shoulder beyond it when that
# shoulder sits < PROX_S from the apex OR the saddle is shallow (valley > MERGE_FRAC of the shoulder);
# stop at the baseline floor or a dead gap. Proximity is apex-relative, so it can't chain to a distant peak.
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
  # baseline & noise from the FULL window (NA/below-detection as 0); detected-points-only inflates noise.
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
  # RP: weight AUC toward the window centre so the on-RT peak beats a larger-area interferent. HILIC: area-select.
  auc <- vapply(peaks, function(p) p$auc, numeric(1))
  score <- if (hilic) auc else
    auc * exp(-0.5 * ((vapply(peaks, function(p) p$apex_rt, numeric(1)) - (acc_lo + acc_hi) / 2) / RT_SIGMA)^2)
  peaks[[which.max(score)]]
}

# Chromatography-conditional accept gate. Shared: point count, intensity floor, SNR, FWHM band, slope.
# HILIC adds return-to-baseline (prominence) + no much-taller neighbour. RP drops those two (they killed
# sharp truncated-flank peaks beside the sulfone) but still makes a 3-point peak clear its own baseline.
accept_peak <- function(p, hilic) {
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
# MsQuality measures the peak the picker chose, handed its boundaries, so width, FWHM,
# prominence and shape are mzQC controlled-vocabulary quantities rather than local
# re-implementations. gaussian_similarity needs >= 5 points and is NA on the sparsest
# peaks, so it is recorded but never gated on. Apex and SNR stay local: both belong to the
# peak the picker selected, and MsQuality cannot be told which peak that is.
add_msquality <- function(chr, picks) {
  n <- length(picks)
  pb <- matrix(NA_real_, n, 2,
               dimnames = list(NULL, c("left_boundary", "right_boundary")))
  for (i in seq_len(n))
    if (!is.null(picks[[i]])) pb[i, ] <- c(picks[[i]]$lb, picks[[i]]$rb)
  if (!any(is.finite(pb[, 1]))) return(picks)
  fwhm  <- MsQuality::xicFwhm(chr, peakBoundary = pb)
  width <- MsQuality::peakWidth(chr, peakBoundary = pb)
  prom  <- MsQuality::peakProminence(chr, peakBoundary = pb)
  gauss <- MsQuality::gaussianSimilarity(chr, peakBoundary = pb)[, "gaussian_similarity"]

  # maxIntensity() takes no peakBoundary, so scope it by handing it a Chromatograms holding
  # only the picked regions. Bit-identical to the local apex on every reverse-phase peak
  # tested; it can differ on a merged HILIC cluster, where the integrated region contains a
  # shoulder taller than the apex the picker chose, so the local value wins there.
  pd  <- Chromatograms::peaksData(chr)
  sub <- lapply(seq_len(n), function(i) {
    if (is.null(picks[[i]])) return(data.frame(rtime = numeric(), intensity = numeric()))
    d <- as.data.frame(pd[[i]])
    d[d$rtime >= picks[[i]]$lb & d$rtime <= picks[[i]]$rb, c("rtime", "intensity")]
  })
  pkchr <- Chromatograms::Chromatograms(
    Chromatograms::ChromBackendMemory(),
    chromData = data.frame(msLevel = 1L, mz = NA_real_,
                           dataOrigin = paste0("p", seq_len(n))),
    peaksData = sub)
  apex <- MsQuality::maxIntensity(pkchr)

  for (i in seq_len(n)) {
    if (is.null(picks[[i]])) next
    # keep the local value where MsQuality cannot return one, so the gate never sees NA
    if (is.finite(fwhm[i]))  picks[[i]]$fwhm    <- fwhm[i]
    if (is.finite(width[i])) picks[[i]]$pkwidth <- width[i]
    # the apex must stay the SELECTED peak's, so adopt MsQuality's only when they agree
    if (is.finite(apex[i]) && isTRUE(all.equal(apex[i], picks[[i]]$apex)))
      picks[[i]]$apex <- apex[i]
    picks[[i]]$msq_prominence <- prom[i]
    picks[[i]]$gauss          <- gauss[i]
  }
  picks
}

snap_one <- function(rt, it, center) {
  ok <- is.finite(rt) & is.finite(it) & it > 0; rt <- rt[ok]; it <- it[ok]
  if (length(rt) < 3) return(NA_real_)
  o <- order(rt); rt <- rt[o]; it <- it[o]
  lm <- MsCoreUtils::localMaxima(it, hws = 1L)
  cand <- which(lm & abs(rt - center) <= SNAP_W); if (!length(cand)) return(NA_real_)
  rt[cand[which.max(it[cand])]]
}

# MS1-apex acceptance window: consensus of snapped confirmed apexes +/- APEX_PAD, then extended
# over peaks contiguous with that span. `bare` normalises a filename to the match key; `conf` is the
# confirmed-hit table, `ov` the manual apex overrides, `old_win` the fallback, `hilic` the picker mode.
acceptance_window <- function(ds, role_full, cache_role, bare, old_win, conf, ov,
                              hilic = FALSE) {
  hh <- conf[conf$dataset_id == ds & conf$role == role_full, ]
  if (!nrow(hh)) return(old_win)
  fk <- bare(cache_role$files); ms2 <- as.numeric(hh$rt)
  ovr <- ov[ov$dataset_id == ds & ov$role == role_full, ]
  apex <- vapply(seq_len(nrow(hh)), function(i) {
    j <- match(bare(hh$source_file[i]), fk); if (is.na(j)) return(NA_real_)
    k <- which(abs(ms2[i] - ovr$ms2_rt) <= 2)
    if (length(k)) return(ovr$apex_rt[k[1]])
    snap_one(cache_role$traces[[j]]$rt, cache_role$traces[[j]]$it, ms2[i])
  }, numeric(1))
  # One anchor per file: a file fragmented 15 times must not out-vote one fragmented 5
  # times, since the number of MS2 events is an artefact of DDA. Take the most intense
  # apex rather than the median, because a file's hits can snap to different peaks and a
  # median would land between them, on nothing.
  keep_a <- is.finite(apex)
  if (!any(keep_a)) return(old_win)
  a_int <- vapply(seq_len(nrow(hh)), function(i) {
    if (!keep_a[i]) return(NA_real_)
    j <- match(bare(hh$source_file[i]), fk); if (is.na(j)) return(NA_real_)
    tr <- cache_role$traces[[j]]
    tr$it[which.min(abs(tr$rt - apex[i]))]
  }, numeric(1))
  a <- vapply(split(seq_len(nrow(hh))[keep_a], bare(hh$source_file)[keep_a]),
              function(ix) apex[ix][which.max(a_int[ix])], numeric(1))
  a <- a[is.finite(a)]; if (!length(a)) return(old_win)
  # Consensus is the best-supported cluster, not the median: anchors are not always
  # unimodal, and a median can fall in the gap between two clusters. Where two clusters are
  # equally supported the deposit does not tell us which retention time is the compound, so
  # both are spanned and the window is flagged uncertain rather than resolved by tie-break.
  nbrs <- vapply(a, function(z) sum(abs(a - z) <= CLUST_TOL), integer(1))
  top  <- a[nbrs == max(nbrs)]
  tied <- diff(range(top)) > CLUST_TOL
  inl  <- if (tied) a else a[abs(a - stats::median(top)) <= CLUST_TOL]
  if (!length(inl)) return(old_win)
  # One pad, never two. The old rule doubled it for a lone anchor, which widens the window exactly
  # when the evidence is thinnest: that is how MSV000082433 admitted a feature 6.9 s off the
  # confirmed RT (37 detections -> 6 once the pad is single).
  w   <- c(min(inl) - APEX_PAD,     max(inl) + APEX_PAD)
  cap <- c(min(inl) - 2 * APEX_PAD, max(inl) + 2 * APEX_PAD)
  # A fixed edge is only fair when it lands in a gap, so extend over peaks that are CONTIGUOUS with
  # what is already accepted (within COHERE_S). That rescues peaks a hard cut would clip 0.4 s from
  # their neighbour, without ever reaching a separate cluster: growth stops at the first real gap
  # and is capped at 2 * APEX_PAD either way.
  cand <- vapply(cache_role$traces, function(t) {
    p <- pick(t$rt, t$it, cap[1], cap[2], hilic)
    if (is.null(p) || !accept_peak(p, hilic)) NA_real_ else p$apex_rt
  }, numeric(1))
  cand <- cand[is.finite(cand)]
  for (i in seq_len(20)) {
    lo <- cand[cand < w[1] & cand >= w[1] - COHERE_S]
    hi <- cand[cand > w[2] & cand <= w[2] + COHERE_S]
    if (!length(lo) && !length(hi)) break
    if (length(lo)) w[1] <- max(min(lo), cap[1])
    if (length(hi)) w[2] <- min(max(hi), cap[2])
  }
  # Record WHY the window is what it is, so a wide one is never mistaken for a confident one.
  attr(w, "n_anchor_files") <- length(a)
  attr(w, "uncertain")      <- tied
  w
}

# Manual reviewer drop-list, keyed on norm_file_key (stable across re-extraction/re-ordering). `drop_all` = ms1_drop_list.csv.
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
          # Pass assayName so the right assay loads: MsBackendMetaboLights (>= 1.7.4, gabri) bakes the
          # assay index into the cache filename, so same-basename pos/neg files no longer collide.
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

# Pick the assay with the most confirmed parent+metabolite features, then derive m/z and RT windows
# from the dominant adduct per role. Returns a list ready for extract_eics().
# rt_pad_s = 5 is load-bearing: 5-OH-omeprazole (CYP2C19) and omeprazole sulfone (CYP3A4) are exact
# isomers at m/z 362.1169, so only RT separates them. In MSV000084008 the 5-OH span is 154.9-158.3 s
# and the sulfone sits at ~181 s; a wide pad lets the sulfone in and it gets integrated instead.
# Validated vs PharmKinetics_CYP2C19: pad 3-20 s r=+0.55..+0.62 (p<0.005), pad 30 s r=-0.01 — do not widen.
derive_extraction_config <- function(dataset_id, conf_final,
                                     rt_pad_s = 5, bio_files = NULL) {
  pc <- conf_final[conf_final$dataset_id == dataset_id, ]
  if (!nrow(pc)) stop("No confirmed features for ", dataset_id)

  # Membership rule: the RT anchor must come only from confirmed hits on files we actually extract.
  # A MASST/SIRIUS hit can sit on a file outside this deposit's set (e.g. 082493's confirmations live
  # on 084008 blood files while 082493 contributes feces/skin/urine); anchoring there imports a foreign
  # matrix's RT. If filtering leaves no parent or metabolite anchor, stop loudly rather than fall back.
  if (!is.null(bio_files)) {
    # norm_file_key(), never a local copy of the rule: confirmation writes back names like
    # ST002044_rawdata_NHF_45.mzML and MTBLS11656_1_MTBLS11656_B-112.mzXML, where the prefix
    # carries an assay segment and can be doubled. A one-prefix strip matched neither, and
    # dropped both deposits at this check.
    in_set <- norm_file_key(pc$source_file, bio_files) %in% norm_file_key(bio_files, bio_files)
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

  # Anchor on the best tier available, per role. This restriction belongs HERE and nowhere else:
  # the RT window should be set by the most confident identification, but a deposit's confirmed
  # feature set must not be thinned just because one of its features reached a higher tier.
  # Applied after the membership and assay filters, so the tier is the best among confirmations
  # that are actually usable for this extraction rather than the best that ever existed.
  if ("confidence" %in% names(ac)) {
    ord <- c(strict = 1, standard = 2, liberal = 3, not_confirmed = 4)
    ac  <- do.call(rbind, lapply(split(ac, ac$role), function(z) {
      r <- ord[as.character(z$confidence)]
      z[!is.na(r) & r == min(r, na.rm = TRUE), ]
    }))
  }

  ad_p  <- dominant_adduct(ac[ac$role == "parent",     ])
  ad_m  <- dominant_adduct(ac[ac$role == "metabolite", ])

  mz_p <- median(ac$mz[ac$role == "parent"     & ac$adduct_inferred == ad_p],
                 na.rm = TRUE)
  mz_m <- median(ac$mz[ac$role == "metabolite" & ac$adduct_inferred == ad_m],
                 na.rm = TRUE)
  # Per-dataset RT window: span this dataset's confirmed RTs and pad by rt_pad_s. The span already
  # absorbs cross-file drift, so the pad only covers files with no confirmed hit — keep it small (isomer note above).
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

# Padding kept in chr_all beyond the detection windows. Must stay >= the widest re-extraction any
# consumer does: dataset-report.qmd inspects "neither" files at +/- 60 s, so 75 s leaves headroom.
CHR_ALL_PAD_S <- 75

extract_eics <- function(ms1, mz_par, mz_met,
                         rt_par_win, rt_met_win, ppm = 20) {
  # Restrict spectra to the RT/mz neighbourhood before building chr_all. Verified on MSV000084008 as
  # AUC-identical, with chr_all much smaller and extraction ~16% faster.
  rt_lo <- min(rt_par_win[1], rt_met_win[1]) - CHR_ALL_PAD_S
  rt_hi <- max(rt_par_win[2], rt_met_win[2]) + CHR_ALL_PAD_S
  mz_lo <- min(mz_par, mz_met) * (1 - 50e-6)
  mz_hi <- max(mz_par, mz_met) * (1 + 50e-6)
  ms1   <- filterMzRange(filterRt(ms1, c(rt_lo, rt_hi)), c(mz_lo, mz_hi))

  all_orig <- unique(dataOrigin(ms1))
  peak_tbl <- rbind(
    # Extract WIDE_EXTRA s past each window so pick()'s valley walk has data; only in-window apexes are accepted.
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
  # chromExtract yields a lazily Spectra-backed object, so every downstream peaksData()/intensity()
  # re-materialises it. setBackend() pulls the targeted EICs into RAM once (count_peaks/pick/plots then
  # read from memory). chr_all stays Spectra-backed — only reports/panels re-chromExtract it, and it stays small.
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

# Per file: build the acceptance window, pick the dominant peak, gate it (accept_peak), apply the drop-list.
# metrics_* keep the auc/fwhm/prominence/pkwidth/snr schema PLUS a `detected` flag (+ apex RT, npts) so
# build_sample_table and 05 read one authoritative detection. Wide per-file EICs kept as chr_*_narrow for reports.
# gauss = MsQuality gaussian_similarity, recorded but never gated on (NA below 5 points).
# msq_prom = MsQuality peakProminence, a peak-to-baseline ratio -- not the same quantity as
# `prominence`, an absolute rise above the higher boundary, and not its threshold.
ROLE_COLS <- c("auc", "fwhm", "prominence", "pkwidth", "snr", "apex_rt", "npts",
               "gauss", "msq_prom", "detected")

process_eics <- function(chr_par, chr_met, chr_all, dataset_id, cfg, conf, ov, drop_all,
                         is_hilic, files_base) {
  RT_ANCHOR <- list()
  peaks_par <- count_peaks(chr_par)
  peaks_met <- count_peaks(chr_met)
  bare <- function(x) norm_file_key(to_bare_name(basename(as.character(x)), files_base))

  role_metrics <- function(chr, role, role_full, old_win) {
    files  <- basename(dataOrigin(chr))
    traces <- lapply(Chromatograms::peaksData(chr),
                     function(pd) list(rt = pd[, "rtime"], it = pd[, "intensity"]))
    win    <- acceptance_window(dataset_id, role_full, list(files = files, traces = traces),
                      bare, old_win, conf, ov, is_hilic)
    RT_ANCHOR[[role]] <<- list(win = as.numeric(win),
                               n_anchor_files = attr(win, "n_anchor_files"),
                               uncertain = isTRUE(attr(win, "uncertain")))
    picks  <- lapply(traces, function(tr) pick(tr$rt, tr$it, win[1], win[2], is_hilic))
    picks  <- add_msquality(chr, picks)
    keep   <- which(!vapply(picks, is.null, logical(1)))
    cols   <- paste0(role, "_", ROLE_COLS)
    if (!length(keep)) {
      m <- data.frame(file = character(0), stringsAsFactors = FALSE)
      for (cc in cols) m[[cc]] <- if (endsWith(cc, "detected")) logical(0) else numeric(0)
      return(m)
    }
    pk <- picks[keep]; fl <- files[keep]
    det <- vapply(seq_along(pk), function(i)
      accept_peak(pk[[i]], is_hilic) &&
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
      gauss      = vapply(pk, function(p) if (is.null(p$gauss)) NA_real_ else p$gauss, 0),
      msq_prom   = vapply(pk, function(p)
                     if (is.null(p$msq_prominence)) NA_real_ else p$msq_prominence, 0),
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
       metrics_par = metrics_par, metrics_met = metrics_met,
       rt_anchor = RT_ANCHOR)
}

# ---- 5. sample-table assembly --------------------------------------------

# NOTE: to_bare_name() lives in R/metadata.R next to norm_file_key() (which calls it); several scripts
# source metadata.R alone, so it can't depend on ms1.R.

# Reject a peak that is at the WRONG COMPOUND'S retention time.
#
# Both target masses carry a known isobar. The metabolite shares 362.1169 with
# omeprazole sulfone, and -- less obviously -- the parent shares 346.1220 with
# hydroxyomeprazole sulfide, which the source study annotates in these very samples
# and which yields the same product ion (m/z 198.0571) as omeprazole, so it also
# passes spectral confirmation. MS1 cannot separate either pair; only retention can.
#
# The two compounds are chromatographically resolved by construction, so a peak in
# one channel sitting nearer the OTHER channel's confirmed retention time is not the
# compound that channel is for. Two situations are handled:
#
# It fires only where both channels' anchors are resolved and far enough apart to
# tell one from the other. Where either is `uncertain` -- its confirmed spectra
# split into two retention clusters -- there is no trustworthy centre to compare
# against, so nothing is rejected and the ambiguity is reported instead.
#
# Applied post hoc on stored apex times, so it costs no re-extraction.
resolve_isobars <- function(metrics_par, metrics_met, conf, dataset_id,
                            rt_anchor = NULL) {
  centre <- function(role_full) {
    q <- as.numeric(
      conf$rt[conf$dataset_id == dataset_id & conf$role == role_full])
    q <- q[is.finite(q)]
    if (length(q)) stats::median(q) else NA_real_
  }
  cp <- centre("parent"); cm <- centre("metabolite")
  # The rule needs BOTH anchors to be real. Where either channel's confirmed
  # spectra split into two retention clusters, its median is a point between two
  # peaks rather than a retention time, and comparing against it rejects peaks
  # that sit exactly where the compound elutes -- which is what it did: it
  # discarded every confirmed-file detection in MSV000088255 and ST002722,
  # including peaks within 4 s of their own anchor. An unresolved anchor is a
  # reason to withhold judgement, not to reject.
  resolved <- !isTRUE(rt_anchor$par$uncertain) && !isTRUE(rt_anchor$met$uncertain) &&
    is.finite(cp) && is.finite(cm) && abs(cp - cm) >= CHAN_SEP
  drop <- function(m, role) {
    col <- paste0(role, "_apex_rt"); det <- paste0(role, "_detected")
    if (!resolved || !all(c(col, det) %in% names(m))) return(m)
    own <- if (role == "par") cp else cm
    oth <- if (role == "par") cm else cp
    bad <- m[[det]] %in% TRUE & is.finite(m[[col]]) &
      abs(m[[col]] - own) >= abs(m[[col]] - oth)
    m[[det]][bad] <- FALSE
    m
  }
  list(par = drop(metrics_par, "par"), met = drop(metrics_met, "met"))
}

# ---- ablations: re-run the extraction with one choice changed ---------------
# `role_metrics()` above is the pipeline; these re-run it from the cached EICs
# with a single parameter moved, so a notebook can ask what any one decision is
# worth. Same picker, same gate, same drop-list -- only the named argument
# differs, and ablation_variant() with no arguments must reproduce the cached
# detection flags exactly.

ROLE_FULL <- c(par = "parent", met = "metabolite")

ablation_setup <- function(ds, bio_files) {
  r  <- readRDS(sprintf("artifacts/per_dataset/%s.rds", ds))
  fb <- basename(as.character(bio_files))
  ch <- list(par = r$chr_par_narrow, met = r$chr_met_narrow)
  list(ds = ds, cfg = r$cfg, ch = ch,
       metrics = list(par = r$metrics_par, met = r$metrics_met),
       tr = lapply(ch, function(c0) lapply(Chromatograms::peaksData(c0),
              function(pd) list(rt = pd[, "rtime"], it = pd[, "intensity"]))),
       fl = lapply(ch, function(c0) basename(dataOrigin(c0))),
       bare = function(x) norm_file_key(to_bare_name(basename(as.character(x)), fb)))
}

# The consensus the pipeline does NOT use: every confirmed hit votes, so a file
# fragmented fifteen times outweighs one fragmented five.
ablation_window_per_hit <- function(e, role, conf, old_win) {
  hh <- conf[conf$dataset_id == e$ds & conf$role == ROLE_FULL[[role]], ]
  a  <- suppressWarnings(as.numeric(hh$rt)); a <- a[is.finite(a)]
  if (!length(a)) return(old_win)
  nbrs <- vapply(a, function(z) sum(abs(a - z) <= CLUST_TOL), integer(1))
  top  <- a[nbrs == max(nbrs)]
  inl  <- a[abs(a - stats::median(top)) <= CLUST_TOL]
  w <- c(min(inl) - APEX_PAD, max(inl) + APEX_PAD)
  attr(w, "uncertain") <- diff(range(top)) > CLUST_TOL
  w
}

# Stops where role_metrics() stops, before the isobar rule; ablation_isobars()
# is that step, kept separate so the baseline can be checked against the cache.
ablation_variant <- function(e, conf, floor_val = MIN_INT, per_hit = FALSE,
                             ov, drop_all) {
  anchor <- list(); mets <- list()
  for (role in c("par", "met")) {
    old_win <- e$cfg[[paste0("rt_", role, "_win")]]
    win <- if (per_hit) ablation_window_per_hit(e, role, conf, old_win) else
      acceptance_window(e$ds, ROLE_FULL[[role]],
                        list(files = e$fl[[role]], traces = e$tr[[role]]),
                        e$bare, old_win, conf, ov, FALSE)
    anchor[[role]] <- list(win = as.numeric(win),
                           uncertain = isTRUE(attr(win, "uncertain")))
    picks <- lapply(e$tr[[role]], function(t) pick(t$rt, t$it, win[1], win[2], FALSE))
    picks <- add_msquality(e$ch[[role]], picks)
    keep  <- which(!vapply(picks, is.null, logical(1)))
    pk <- picks[keep]; fl <- e$fl[[role]][keep]
    old_floor <- MIN_INT; MIN_INT <<- floor_val
    det <- vapply(seq_along(pk), function(i) accept_peak(pk[[i]], FALSE) &&
                    !is_dropped(e$ds, ROLE_FULL[[role]], fl[i], drop_all), logical(1))
    MIN_INT <<- old_floor
    mets[[role]] <- stats::setNames(
      data.frame(file = fl,
                 auc = vapply(pk, function(p) p$auc, numeric(1)),
                 apex_rt = vapply(pk, function(p) p$apex_rt, numeric(1)),
                 detected = det, stringsAsFactors = FALSE),
      c("file", paste0(role, c("_auc", "_apex_rt", "_detected"))))
  }
  list(ds = e$ds, win = anchor, par = mets$par, met = mets$met)
}

ablation_isobars <- function(v, conf) {
  mr <- resolve_isobars(v$par, v$met, conf, v$ds, v$win)
  modifyList(v, list(par = as.data.frame(mr$par), met = as.data.frame(mr$met)))
}

build_sample_table <- function(metrics_par, metrics_met, files_base, meta,
                               dedup = TRUE) {
  # dedup = TRUE -> one row per subject (default; 05 clinical groups, avoids pseudoreplication).
  # dedup = FALSE -> one row per file (subject_id retained) for genotype/longitudinal/specimen analyses.
  # Merge on norm_file_key, never the bare name: some cache names carry the accession twice
  # (MTBLS1866_1_MTBLS1866_DA44_p.mzML), which no bare-name match resolves. That silently
  # demoted 19 detected files to "neither" before it was caught.
  metrics_par$.k <- norm_file_key(to_bare_name(metrics_par$file, files_base))
  metrics_met$.k <- norm_file_key(to_bare_name(metrics_met$file, files_base))
  metrics_par$file <- NULL; metrics_met$file <- NULL
  sample_tbl <- merge(metrics_par, metrics_met, by = ".k", all = TRUE)
  sample_tbl <- merge(data.frame(file = files_base, .k = norm_file_key(files_base),
                                 stringsAsFactors = FALSE),
                      sample_tbl, by = ".k", all.x = TRUE)
  sample_tbl$.k <- NULL

  # Join metadata on the NORMALISED key, not the raw name: a source may spell the extension differently
  # (084008's submitter TSV says .mzXML, the deposit .mzML), so a literal merge matched 0/461 rows and
  # left subject_id/body_site/timepoint NA for the whole cohort. norm_file_key drops extension+accession (439 match).
  sample_tbl$.k <- norm_file_key(sample_tbl$file)
  meta          <- meta[!duplicated(norm_file_key(meta$file)), , drop = FALSE]
  meta$.k       <- norm_file_key(meta$file)
  meta$file     <- NULL                     # keep the deposit's spelling
  sample_tbl    <- merge(sample_tbl, meta, by = ".k", all.x = TRUE)
  sample_tbl$.k <- NULL

  # Detection is decided once in the picker (accept_peak + drop-list); metrics_* carry par_detected/
  # met_detected. Files with no picked peak merge in as NA -> FALSE. Gate documented in 03-all-datasets.qmd.
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

  # Dedup to one file per subject by highest combined AUC. Fall back to sample_name then file when
  # subject_id is NA, else all-NA rows collapse into one group.
  # Take the key columns defensively: a metadata source can fail (the Workbench fetch is a
  # network call) and hand back a frame with only `file`, in which case subject_id/sample_name
  # are absent rather than NA and coalesce() dies on a zero-length vector.
  col_or_na <- function(nm) if (nm %in% names(sample_tbl)) as.character(sample_tbl[[nm]])
                            else rep(NA_character_, nrow(sample_tbl))
  sample_tbl$.dedup_key <- dplyr::coalesce(
    col_or_na("subject_id"),
    col_or_na("sample_name"),
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

# Cross-dataset entry point. Resolves deposit_id/assay/mz/RT windows from the artifacts; caller passes dataset_id + three artifacts.
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
  # For MetaboLights, pass the assay so same-basename pos/neg modes are disambiguated at fetch time.
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
    # Per role: the window actually used, how many FILES anchored it, and whether the
    # anchors were split between equally-supported retention times. A wide window with
    # uncertain = TRUE is an unresolved anchor, not a confident wide peak.
    rt_anchor      = proc$rt_anchor,
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

# Shared y-axis limit across all EICs; 5% headroom keeps the apex off the border.
eic_ylim <- function(chr) {
  ints <- unlist(intensity(chr), use.names = FALSE)
  c(0, max(ints, na.rm = TRUE) * 1.05)
}

# ===== peak detection + metrics (was peaks.R) =====

# Count peaks per EIC via a local-maxima detector. hws: half-window in scans (max within ±hws
# neighbours). pct_max: min height as a fraction of the trace max (filters sub-noise wiggles).
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

# Flag files whose FWHM deviates > 3 × MAD from the dataset median. NA inputs (no detected peak)
# count as non-outliers so they don't propagate into downstream sum() calls.
flag_fwhm <- function(df, col) {
  v   <- df[[col]]
  med <- median(v, na.rm = TRUE)
  mad <- stats::mad(v, na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) return(rep(FALSE, length(v)))
  flagged <- abs(v - med) > 3 * mad
  flagged[is.na(flagged)] <- FALSE
  flagged
}

# Report whether SIRIUS-confirmed files produced a peak. They are positive controls: a miss is a
# settings problem (ppm/RT/adduct), not biology. Extension-agnostic match (.mzML vs .mzXML re-deposits).
check_confirmed <- function(files, n_peaks, conf_files, role) {
  # Canonical rule, not a private copy: a local one-prefix regex is what silently lost 21 of
  # MTBLS1866's 23 anchors elsewhere. Both sides are cache-derived here, so they agree today,
  # but the point is not to keep a second implementation that can drift.
  idx <- which(norm_file_key(files) %in% norm_file_key(conf_files))
  cat(role, "- confirmed files loaded:", length(idx),
      "| with >=1 peak:", sum(n_peaks[idx] > 0), "\n")

  missed <- files[idx[n_peaks[idx] == 0]]
  if (length(missed))
    cat("  WARNING missed:", paste(missed, collapse = ", "), "\n")
}

# Description, not inference: detection tier by group and the per-group ratio summary.
# The per-dataset Fisher and Kruskal tests that used to live here tested files rather than
# subjects, so their p-values estimated nothing. Contrast testing belongs in `05`, under
# the analysis-unit rule.
describe_by_group <- function(sample_tbl, group_col = "health_status") {
  if (!group_col %in% names(sample_tbl)) {
    warning(group_col, " not in sample_tbl - skipping group description")
    return(NULL)
  }
  g <- sample_tbl[[group_col]]
  both <- sample_tbl[sample_tbl$detection_tier == "both" &
                     !sample_tbl$par_fwhm_outlier %in% TRUE &
                     !sample_tbl$met_fwhm_outlier %in% TRUE, ]
  per_group_ratio <- if (nrow(both) > 0)
    dplyr::group_by(both, .data[[group_col]]) |>
    dplyr::summarise(n          = dplyr::n(),
                     median_l2r = median(log2_ratio, na.rm = TRUE),
                     iqr_l2r    = stats::IQR(log2_ratio, na.rm = TRUE),
                     .groups    = "drop") else NULL
  list(group_col       = group_col,
       detection_table = table(detection_tier = sample_tbl$detection_tier,
                               group = g, useNA = "ifany"),
       per_group_ratio = per_group_ratio)
}

# Canonical fields usable as stratifiers: categorical, present, >= min_levels values each backed by
# >= min_per_level samples. `exclude` drops fields handled elsewhere; id/name/file/age never qualify.
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

# Plot FWHM-outlier EICs alongside a reference (the file nearest the median FWHM). A genuine outlier
# looks deformed vs the reference; if it looks similar the flag is a false alarm and the file stays in.
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
