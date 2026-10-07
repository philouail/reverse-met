# R/ms1.R: MS1 pipeline, from step 4's ion, m/z and RT range to a per-file sample table:
# load_spectra_ms1 -> extract_eics -> process_eics -> build_sample_table (run_ms1_for_dataset).

# ==========================================================================
# ---- 0. peak-picker core -------------------------------------------------
# One file's (rt, it) + acceptance window -> one integrated peak (or NULL), the accept gate
# and the apex-window builder. Thresholds are explained in 05-extraction.qmd.

WIDE_EXTRA <- 50      # s of EIC read past step 4's range, for the valley walk
HWS        <- 1L      # local-maxima half-window (scans)
GAP_FACTOR <- 4       # stop the walk at a gap > this many median scan spacings
GAP_MIN_S  <- 15      # ... and never below this many s
FWHM_BAND_HILIC <- c(1.3, Inf)  # HILIC FWHM band (s)
FWHM_BAND_RP    <- c(1, 30)     # RP FWHM band (s)
MERGE_FRAC <- 0.65    # merge a saddle if its valley > this frac of the lower top
PROX_S     <- 4       # ... and always merge a maximum < this many s from the apex
NEAR_S     <- 15      # taller-neighbour search radius (s)
NBR_MAX    <- 1.3     # ... reject if a neighbour is > this x the apex
PROM_MIN   <- 0.25    # apex rise above its higher boundary >= this frac of apex
SNR_MIN_HILIC <- 3    # HILIC SNR floor (apex above local baseline, noise units)
SNR_MIN_RP    <- 12   # RP SNR floor
MIN_INT       <- 500  # apex-intensity floor, both modes
MAXW_FWHM  <- 3       # cap each boundary at this many FWHM from the apex
TRUNC_DIP  <- 0.20    # truncated peak: allowed dip in rise/fall, frac of apex
APEX_PAD   <- 5       # s pad of acceptance window; also tail-rejection tolerance
RT_SIGMA   <- 4       # RP: gaussian sd (s) of AUC weight around window centre
SNAP_W     <- 20      # max s from an MS2 rt to its MS1 apex, else no anchor
CLUST_TOL  <- 10      # s: anchors kept around the best-supported cluster
COHERE_S   <- 1       # s: extend the window over peaks this close to its edge
CHAN_SEP   <- 6       # s: closer anchors make resolve_isobars() stand down

# Walk one flank from the apex to its boundary index: descend to a valley, merge the shoulder
# beyond it if within PROX_S of the apex or the saddle is shallow (MERGE_FRAC); stop at floor or gap.
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
  # baseline & noise from the FULL trace (NA/below-detection as 0); detected-points-only inflates noise.
  it0 <- it; it0[!is.finite(it0) | it0 < 0] <- 0
  baseline <- as.numeric(stats::quantile(it0, 0.10))
  bpts  <- it0[it0 <= stats::median(it0)]
  noise0 <- if (length(bpts) >= 3) stats::mad(bpts, constant = 1.4826) else NA_real_
  noise_at <- function(a) if (!is.finite(noise0) || noise0 <= 0) 0.01 * max(a - baseline, 1) else noise0
  rt_full <- rt                     # every scan, empty ones included (truncation test, auc0)
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
    # Truncated peak: all kept points above half height and the scans past both edges empty (the ion left
    # the m/z window). Width runs to those scans at 0; needs contiguous unimodal scans (dip <= TRUNC_DIP).
    truncated <- FALSE
    if (!is.finite(fwhm) && length(ab) && min(ab) == li && max(ab) == ri) {
      jl <- match(rt[li], rt_full) - 1L; jr <- match(rt[ri], rt_full) + 1L
      inside   <- which(rt_full >= rt[li] & rt_full <= rt[ri])
      # the deepest dip: a drop below the highest scan so far on the way up, or a
      # rise above the lowest so far on the way down
      rise <- it[li:ai]; fall <- it[ai:ri]
      dip  <- max(cummax(rise) - rise, fall - cummin(fall))
      unimodal <- dip <= TRUNC_DIP * apex
      if (!is.na(jl) && !is.na(jr) && jl >= 1L && jr <= length(rt_full) &&
          it0[jl] == 0 && it0[jr] == 0 && all(it0[inside] > 0) && unimodal) {
        fwhm <- approx(c(it[ri], 0), c(rt[ri], rt_full[jr]), half)$y -
                approx(c(0, it[li]), c(rt_full[jl], rt[li]), half)$y
        truncated <- TRUE
      }
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
    # Area: trapezoids above a straight baseline from yl at the left boundary to yr at the right.
    line_at <- function(yl, yr) if (rt[ri] > rt[li]) yl + (yr - yl) * (rt[idx] - rt[li]) / (rt[ri] - rt[li])
                                else rep(min(yl, yr), npts)
    area_on <- function(bl) {
      ic <- pmax(it[idx] - bl, 0)
      sum(diff(rt[idx]) * (head(ic, -1) + tail(ic, -1)) / 2)
    }
    bl  <- line_at(it[li], it[ri])   # the line between the two boundary scans
    auc <- area_on(bl)
    # auc0: the same area with the line from 0 at a boundary whose next scan is empty; only step 6's
    # baseline check uses it. Boundaries, SNR, gate and selection stay on `bl` and `auc`.
    jl0 <- match(rt[li], rt_full) - 1L; jr0 <- match(rt[ri], rt_full) + 1L
    yl0 <- if (jl0 >= 1L && it0[jl0] == 0) 0 else it[li]
    yr0 <- if (jr0 <= length(rt_full) && it0[jr0] == 0) 0 else it[ri]
    auc0 <- area_on(line_at(yl0, yr0))
    apex_base <- bl[ai - li + 1L]
    nb <- lm & abs(rt - a_rt) <= NEAR_S & !(rt >= rt[li] & rt <= rt[ri])
    nbr_ratio <- if (any(nb)) max(it[nb]) / apex else 0
    list(apex_rt = a_rt, apex = apex, lb = rt[li], rb = rt[ri], idx = idx,
         lb_int = it[li], rb_int = it[ri], npts = npts, auc = auc, auc0 = auc0, fwhm = fwhm,
         truncated = truncated,
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

# Accept gate: gate_reason() names the first rule a pick fails, or "accepted". Both modes test intensity,
# SNR and FWHM band; HILIC adds return to baseline and no taller neighbour; RP tests 3-point peaks' baseline.
gate_reason <- function(p, hilic) {
  if (is.null(p)) return("no peak in window")
  if (!is.finite(p$auc)) return("no area")
  fw      <- if (hilic) FWHM_BAND_HILIC else FWHM_BAND_RP
  snr_min <- if (hilic) SNR_MIN_HILIC  else SNR_MIN_RP
  if (p$apex < MIN_INT)                      return("below intensity floor")
  if (!is.finite(p$snr) || p$snr < snr_min)  return("SNR below floor")
  if (!is.finite(p$fwhm))                    return("FWHM not measurable")
  if (p$fwhm < fw[1] || p$fwhm > fw[2])      return("FWHM outside band")
  baseline_ok <- p$prominence >= PROM_MIN * p$apex
  if (hilic) {
    if (!baseline_ok)          return("no return to baseline")
    if (p$nbr_ratio > NBR_MAX) return("taller neighbour")
  } else if (p$npts < 4L && !baseline_ok) return("3-point peak off baseline")
  "accepted"
}
accept_peak <- function(p, hilic) identical(gate_reason(p, hilic), "accepted")

# Replace the picker's FWHM and width with MsQuality's on the picked boundaries; add its prominence and
# gaussian similarity (NA under 5 points, never gated). Apex and SNR stay local, from the selected peak.
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

  # maxIntensity() takes no peakBoundary, so hand it Chromatograms holding only the picked regions.
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

# MS1 apex rt for one confirmed MS2 rt (`center`): walk both flanks from the nearest scan as pick() does
# and take the region's apex, so a spectrum on a tail anchors on its own peak. NA if beyond SNAP_W.
snap_one <- function(rt, it, center) {
  o <- order(rt); rt <- rt[o]; it <- it[o]
  gap_thr <- max(GAP_FACTOR * stats::median(diff(rt)), GAP_MIN_S)
  it0 <- it; it0[!is.finite(it0) | it0 < 0] <- 0
  baseline <- as.numeric(stats::quantile(it0, 0.10))
  bpts   <- it0[it0 <= stats::median(it0)]
  noise0 <- if (length(bpts) >= 3) stats::mad(bpts, constant = 1.4826) else NA_real_
  ok <- is.finite(it) & it > 0; rt <- rt[ok]; it <- it[ok]
  if (length(rt) < 3) return(NA_real_)
  s     <- which.min(abs(rt - center))
  noise <- if (!is.finite(noise0) || noise0 <= 0) 0.01 * max(it[s] - baseline, 1) else noise0
  floor_lv <- baseline + 3 * noise
  li <- walk_bound(it, s, -1L, floor_lv, it[s], rt, gap_thr, PROX_S)
  ri <- walk_bound(it, s, +1L, floor_lv, it[s], rt, gap_thr, PROX_S)
  a  <- (li:ri)[which.max(it[li:ri])]
  if (abs(rt[a] - center) > SNAP_W) return(NA_real_)
  rt[a]
}

# MS1 acceptance window: snapped confirmed apexes +/- APEX_PAD, extended over contiguous peaks. `bare` makes
# the file match key; `conf` confirmed hits, `ov` apex overrides, `old_win` fallback, `hilic` picker mode.
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
  # One anchor per file, its most intense apex: MS2 event counts are a DDA artefact, and a median
  # of a file's hits could fall between two peaks.
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
  # Consensus is the best-supported cluster (a median can fall between two). Equally supported
  # clusters are all spanned and the window flagged uncertain, not resolved by a tie-break.
  nbrs <- vapply(a, function(z) sum(abs(a - z) <= CLUST_TOL), integer(1))
  top  <- a[nbrs == max(nbrs)]
  tied <- diff(range(top)) > CLUST_TOL
  inl  <- if (tied) a else a[abs(a - stats::median(top)) <= CLUST_TOL]
  if (!length(inl)) return(old_win)
  # One pad, never two: doubling it for a lone anchor would widen the window when evidence is thinnest.
  w   <- c(min(inl) - APEX_PAD,     max(inl) + APEX_PAD)
  cap <- c(min(inl) - 2 * APEX_PAD, max(inl) + 2 * APEX_PAD)
  # Extend each edge over accepted peaks within COHERE_S of it, so a hard cut does not clip a peak;
  # growth stops at the first gap and is capped at 2 * APEX_PAD.
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
  attr(w, "anchors")        <- as.numeric(inl)   # the file anchors the window rests on
  w
}

# Manual reviewer drop-list, keyed on norm_file_key (stable across re-extraction/re-ordering). `drop_all` = ms1_drop_list.csv.
is_dropped <- function(dataset, role_full, file, drop_all) {
  if (!nrow(drop_all)) return(FALSE)
  dl <- drop_all[drop_all$dataset == dataset & drop_all$role == role_full, ]
  nrow(dl) > 0 && norm_file_key(file) %in% norm_file_key(dl$file)
}

# ---- 1. spectrum loading (repo dispatch) ---------------------------------

# Every spectrum of `files` of one deposit through its repository's MsBackend: files are synced into
# the BiocFileCache (not downloaded again) and read offline. load_spectra_ms1() keeps the MS1 scans.
load_spectra <- function(deposit_id, files, assay = character()) {
  wanted <- basename(files)
  pat    <- paste0("(", paste(gsub("\\.", "\\\\.", wanted),
                              collapse = "|"), ")$")
  be <- if (startsWith(deposit_id, "MTBLS")) {
          # Pass assayName so the right assay loads: MsBackendMetaboLights (>= 1.7.4, gabri) bakes the
          # assay index into the cache filename, so same-basename pos/neg files do not collide.
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
  Spectra(be)
}

load_spectra_ms1 <- function(deposit_id, files, assay = character())
  filterMsLevel(load_spectra(deposit_id, files, assay), 1L)

# ---- 1b. instrument and m/z tolerance ------------------------------------
# Step 4 reads each hit file's instrument (read_instrument()) and classes each assay (instrument_class());
# step 5 extracts at that class's tolerance, the wider one when the instrument is unknown.

PPM_BY_CLASS <- c(orbitrap = 5, tof = 20)   # step 5's m/z tolerance (ppm) by mass analyser
PPM_DEFAULT  <- 20                          # ... when the instrument cannot be told

# "orbitrap", "tof" or "unknown", from a file header's (or deposit metadata's) instrument fields.
# The model decides: converters report a Q Exactive's analyser as "quadrupole", its first stage.
instrument_class <- function(model, analyzer = NA, manufacturer = NA) {
  x <- tolower(paste(model, analyzer, manufacturer))
  ifelse(grepl("orbitrap|exactive|exploris", x), "orbitrap",
  ifelse(grepl("tof|time-of-flight|maxis|impact|compact|xevo|synapt", x), "tof", "unknown"))
}

# Local path of each raw file (`files`: cache file names) already in the BiocFileCache, else NA. Uses
# bfcinfo(): the packages' own offline lookups only know the deposit synced last.
cached_raw_paths <- function(files) {
  rpath <- BiocFileCache::bfcinfo(BiocFileCache::BiocFileCache(ask = FALSE))$rpath
  p <- rpath[match(basename(files), basename(rpath))]
  ifelse(!is.na(p) & file.exists(p), p, NA_character_)
}

# One raw file's header: the instrument it records and the polarity of a few of its scans
# ("mixed" where those scans switch polarity).
read_instrument <- function(path, n_scans = 10) {
  ms <- mzR::openMSfile(path)
  on.exit(mzR::close(ms))
  ii  <- mzR::instrumentInfo(ms)
  n   <- length(ms)
  pol <- mzR::header(ms, unique(round(seq(1, n, length.out = min(n_scans, n)))))$polarity
  pol <- unique(c(`1` = "positive", `0` = "negative")[as.character(pol)])
  pol <- pol[!is.na(pol)]
  data.frame(manufacturer = ii$manufacturer, model = ii$model, analyzer = ii$analyzer,
             ion_source   = ii$ionisation,
             polarity     = if (length(pol) == 1) pol else if (length(pol)) "mixed" else NA_character_)
}

# ---- 2. derive the ion, m/z and RT range from confirmed features --------

# The ion to extract: the molecular ion ([M+H]+ or [M-H]-) when confirmed, else the most frequent;
# a plain vote can let an in-source fragment win a near-tie.
dominant_adduct <- function(df) {
  tab <- table(df$adduct_inferred[!is.na(df$adduct_inferred)])
  if (!length(tab)) return(NA_character_)
  mol <- intersect(c("[M+H]+", "[M-H]-"), names(tab))
  if (length(mol)) return(mol[which.max(tab[mol])])
  names(tab)[which.max(tab)]
}

# Per role, keep only the confirmations at the best tier reached (strict, then standard, then
# liberal): the ion, m/z and RT range step 4 writes are set by the most confident identifications.
best_tier <- function(ac) {
  if (!"confidence" %in% names(ac)) return(ac)
  ord <- c(strict = 1, standard = 2, liberal = 3, not_confirmed = 4)
  do.call(rbind, lapply(split(ac, ac$role), function(z) {
    r <- ord[as.character(z$confidence)]
    z[!is.na(r) & r == min(r, na.rm = TRUE), ]
  }))
}

# The assay with every target confirmed and the most features, and per target its ion, m/z and padded
# confirmed-RT range (make_cfg()). The range only bounds the trace read; acceptance_window() sets the window.
derive_extraction_config <- function(dataset_id, conf_final,
                                     rt_pad_s = 5, bio_files = NULL,
                                     roles = compound_spec$role) {
  pc <- conf_final[conf_final$dataset_id == dataset_id, ]
  if (!nrow(pc)) stop("No confirmed features for ", dataset_id)

  # With bio_files, only confirmations on those files set the range (a hit on another deposit's file
  # would import a foreign matrix's RT), and a target left with none stops. Step 4 passes no bio_files.
  if (!is.null(bio_files)) {
    # norm_file_key(), not a local rule: confirmation names can carry a doubled, assay-bearing prefix.
    in_set <- norm_file_key(pc$source_file, bio_files) %in% norm_file_key(bio_files, bio_files)
    check_in_set(pc, in_set, dataset_id, roles)
    pc <- pc[in_set, ]
  }

  assay_tab <- pc |>
    dplyr::group_by(assay) |>
    dplyr::summarise(every   = all(roles %in% role),
                     n_total = sum(role %in% roles),
                     .groups = "drop") |>
    dplyr::filter(every) |>
    dplyr::arrange(desc(n_total))
  if (!nrow(assay_tab))
    stop("No assay with ", targets_phrase(length(roles)), " confirmed in ", dataset_id)

  assay <- assay_tab$assay[1]
  ac    <- if (is.na(assay)) pc[is.na(pc$assay), ]
           else              pc[!is.na(pc$assay) & pc$assay == assay, ]

  # Ion, m/z and range from each target's best usable tier (after the membership and assay filters);
  # only here: the MS1 window still anchors on every tier, so the confirmed set is not thinned.
  ac <- best_tier(ac)

  per_role <- lapply(roles, function(r) {
    ad <- dominant_adduct(ac[ac$role == r, ])
    # Per-dataset RT range: span this dataset's confirmed RTs and pad by rt_pad_s. The span already
    # absorbs cross-file drift; the pad only widens the trace read and the fallback window.
    rt <- ac$rt[ac$role == r & ac$adduct_inferred == ad]
    rt <- rt[is.finite(rt)]
    list(ion   = ad,
         mz    = median(ac$mz[ac$role == r & ac$adduct_inferred == ad], na.rm = TRUE),
         win   = c(min(rt) - rt_pad_s, max(rt) + rt_pad_s),
         files = unique(basename(ac$source_file[ac$role == r & !is.na(ac$source_file)])))
  })
  make_cfg(assay, roles, per_role)
}

# The membership rule's stop: each target needs a confirmation on the files step 5 extracts
# (`in_set`, one TRUE/FALSE per row of `pc`).
check_in_set <- function(pc, in_set, dataset_id, roles) {
  n <- vapply(roles, function(r) sum(in_set & pc$role == r), integer(1))
  if (any(n == 0))
    stop("Membership filter leaves ", dataset_id, " with ",
         paste0(target_prefix(roles), "=", n, collapse = " "),
         " in-set confirmations; needs re-confirmation on ",
         "its own extraction set before extraction.")
}

# The cfg list (cached per assay): assay, then per target in `roles` order adduct_<p>, mz_<p>,
# rt_<p>_win and confirmed_files_<role>, from `per_role` (one list of ion, mz, win, files per target).
make_cfg <- function(assay, roles, per_role) {
  p     <- target_prefix(roles)
  field <- function(k) lapply(per_role, function(x) x[[k]])
  c(list(assay = assay),
    stats::setNames(field("ion"),   paste0("adduct_", p)),
    stats::setNames(field("mz"),    paste0("mz_", p)),
    stats::setNames(field("win"),   paste0("rt_", p, "_win")),
    stats::setNames(field("files"), paste0("confirmed_files_", roles)))
}

# Which confirmations set a range of step 4's table (`windows`): those of the row's assay and
# ion, at the best tier the target reached in that assay. One TRUE/FALSE per row of `conf`.
window_setters <- function(conf, windows) {
  conf$row <- seq_len(nrow(conf))
  rows <- unlist(lapply(split(windows, windows$dataset_id), function(w) {
    ac <- best_tier(conf[conf$dataset_id == w$dataset_id[1] & conf$assay %in% w$assay, ])
    ion <- w$ion[match(ac$role, w$role)]
    ac$row[!is.na(ac$adduct_inferred) & ac$adduct_inferred == ion]
  }))
  conf$row %in% rows
}

# Step 4's ion, m/z and RT range for one assay (`windows`) as a cfg. Each target needs a confirmation
# on `bio_files`; if a range-setting one lies outside them, the range is rebuilt without it.
extraction_window <- function(dataset_id, windows, conf_final, bio_files,
                              roles = compound_spec$role) {
  w    <- windows[windows$dataset_id == dataset_id, , drop = FALSE]
  rows <- lapply(roles, function(r) w[w$role == r, , drop = FALSE])
  if (any(vapply(rows, nrow, integer(1)) != 1))
    stop("No ion, m/z and range for ", dataset_id,
         " in artifacts/extraction_windows.csv; render step 4 first.")
  # The assay and the pad are the same on every row of the assay.
  first <- rows[[1]]

  pc <- conf_final[conf_final$dataset_id == dataset_id, ]
  in_set <- norm_file_key(pc$source_file, bio_files) %in% norm_file_key(bio_files, bio_files)
  check_in_set(pc, in_set, dataset_id, roles)

  # Step 4's range rests on confirmations step 5 does not extract: rebuild it without them.
  n_out <- sum(window_setters(pc, w) & !in_set)
  if (n_out > 0) {
    message("  ", dataset_id, ": step 5's range differs from step 4's, because ", n_out,
            " confirmation(s) on files step 5 does not extract were left out of it")
    cfg <- derive_extraction_config(dataset_id, conf_final, rt_pad_s = first$rt_pad_s,
                                    bio_files = bio_files, roles = roles)
    cfg$assay <- as.character(cfg$assay)   # the same type as below
    return(cfg)
  }

  ac <- best_tier(pc[in_set & pc$assay %in% first$assay, ])
  files <- function(role)
    unique(basename(as.character(ac$source_file[ac$role == role & !is.na(ac$source_file)])))
  per_role <- lapply(seq_along(roles), function(i)
    list(ion = rows[[i]]$ion, mz = rows[[i]]$mz, win = c(rows[[i]]$rt_lo, rows[[i]]$rt_hi),
         files = files(roles[i])))
  # A MassIVE assay has no assay name: the column reads back as logical NA.
  make_cfg(as.character(first$assay), roles, per_role)
}

# ---- 3. EIC extraction ---------------------------------------------------

# Padding kept in chr_all beyond step 4's ranges, so a trace can be inspected around the
# window (R/peak_review_panels.R).
CHR_ALL_PAD_S <- 75

# One trace per file and target (`roles`, in compound.csv order), at the m/z and over the RT
# range of `cfg`: chr_all, then chr_<prefix> per target (chr_par, chr_met).
extract_eics <- function(ms1, cfg, ppm = PPM_DEFAULT, roles = compound_spec$role) {
  p   <- target_prefix(roles)
  mz  <- lapply(p, function(x) cfg[[paste0("mz_", x)]])
  win <- lapply(p, function(x) cfg[[paste0("rt_", x, "_win")]])
  # Keep only the RT/m/z neighbourhood of all targets before building chr_all (AUCs unchanged).
  rt_lo <- min(vapply(win, function(w) w[1], numeric(1))) - CHR_ALL_PAD_S
  rt_hi <- max(vapply(win, function(w) w[2], numeric(1))) + CHR_ALL_PAD_S
  mz_lo <- min(unlist(mz)) * (1 - 50e-6)
  mz_hi <- max(unlist(mz)) * (1 + 50e-6)
  ms1   <- filterMzRange(filterRt(ms1, c(rt_lo, rt_hi)), c(mz_lo, mz_hi))

  all_orig <- unique(dataOrigin(ms1))
  # Extract WIDE_EXTRA s past each range so pick()'s valley walk has data; pick() starts only from maxima in the MS1 window.
  peak_tbl <- do.call(rbind, lapply(seq_along(roles), function(i)
    data.frame(role = roles[i], msLevel = 1L, dataOrigin = all_orig,
               rtMin = win[[i]][1] - WIDE_EXTRA, rtMax = win[[i]][2] + WIDE_EXTRA,
               mzMin = mz[[i]] * (1 - ppm * 1e-6),
               mzMax = mz[[i]] * (1 + ppm * 1e-6),
               stringsAsFactors = FALSE)))
  chr_all  <- Chromatograms(ms1)
  # setBackend() pulls the lazily Spectra-backed EICs into RAM once, so later reads are from memory;
  # chr_all stays Spectra-backed (only the review panels re-extract it).
  chr_eics <- setBackend(
    chromExtract(chr_all, peak_tbl, by = c("msLevel", "dataOrigin")),
    ChromBackendMemory())
  c(list(chr_all = chr_all),
    stats::setNames(lapply(roles, function(r) chr_eics[chromData(chr_eics)$role == r]),
                    paste0("chr_", p)))
}

# ---- 3b. raw centroids, for 05's tolerance check ---------------------------
# An EIC keeps no m/z, so the check reads the raw MS1 centroids; the parallel workers'
# backend messages are dropped.

# Per MS1 scan between rt[1] and rt[2], the most intense centroid within `ppm` of `mz`:
# its deviation from `mz` (ppm) and its intensity. A scan with no centroid there is left out.
scan_centroids <- function(ms1, mz, rt, ppm = PPM_DEFAULT) {
  s   <- filterMzRange(filterRt(ms1, rt), mz * (1 + c(-1, 1) * ppm * 1e-6))
  pk  <- suppressMessages(peaksData(s))
  has <- vapply(pk, nrow, integer(1)) > 0
  top <- vapply(pk[has], function(p) p[which.max(p[, "intensity"]), ], numeric(2))
  data.frame(file = basename(dataOrigin(s))[has], rt = rtime(s)[has],
             dev = 1e6 * (top[1, ] - mz) / mz, int = top[2, ])
}

# For the mass-calibration check: ions (centroids >= `min_int` in rt, binned to 0.01 m/z) found in
# >= `min_frac` of the files, with each file's intensity-weighted m/z per ion.
common_ions <- function(ms1, rt, min_int = 1e5, min_frac = 0.9) {
  s  <- filterIntensity(filterRt(ms1, rt), intensity = c(min_int, Inf))
  pk <- suppressMessages(peaksData(s))
  d  <- data.frame(file = rep(basename(dataOrigin(s)), vapply(pk, nrow, integer(1))),
                   mz   = unlist(lapply(pk, function(p) p[, "mz"])),
                   int  = unlist(lapply(pk, function(p) p[, "intensity"])))
  d$ion   <- round(d$mz, 2)
  n_files <- length(unique(basename(dataOrigin(s))))
  d |>
    dplyr::group_by(ion) |>
    dplyr::filter(dplyr::n_distinct(file) >= min_frac * n_files) |>
    dplyr::group_by(file, ion) |>
    dplyr::summarise(mz = sum(mz * int) / sum(int), .groups = "drop") |>
    as.data.frame()
}

# ---- 4. peak pick (per file) + metrics + FWHM flag -----------------------

# Per-file metric columns per target (<prefix>_<col>). auc and auc0 are NA on a truncated peak; gauss,
# msq_prom and auc0 are never gated on; gate = gate_reason(), "manual drop/accept" or "isobar rule".
ROLE_COLS <- c("auc", "auc0", "fwhm", "prominence", "pkwidth", "snr", "apex_rt", "npts",
               "gauss", "msq_prom", "detected", "truncated", "gate")

# `chr` holds each target's traces, named by prefix (per_target(eics, "chr")). Returns, per
# target in that order, peaks_<prefix>, chr_<prefix>_narrow and metrics_<prefix>, then rt_anchor.
process_eics <- function(chr, chr_all, dataset_id, cfg, conf, ov, drop_all,
                         is_hilic, files_base, accept = curation_inputs()$accept) {
  RT_ANCHOR <- list()
  peaks <- lapply(chr, count_peaks)
  bare <- function(x) norm_file_key(to_bare_name(basename(as.character(x)), files_base))

  role_metrics <- function(chr, role, role_full, old_win) {
    files  <- basename(dataOrigin(chr))
    traces <- lapply(Chromatograms::peaksData(chr),
                     function(pd) list(rt = pd[, "rtime"], it = pd[, "intensity"]))
    win    <- acceptance_window(dataset_id, role_full, list(files = files, traces = traces),
                      bare, old_win, conf, ov, is_hilic)
    RT_ANCHOR[[role]] <<- list(win = as.numeric(win),
                               n_anchor_files = attr(win, "n_anchor_files"),
                               uncertain = isTRUE(attr(win, "uncertain")),
                               anchors = attr(win, "anchors"))
    picks  <- lapply(traces, function(tr) pick(tr$rt, tr$it, win[1], win[2], is_hilic))
    # A peak accepted by hand (curation/ms1_accept_list.csv) is picked around the
    # reviewer's apex instead of in the window; the peak criteria still apply.
    acc  <- accept[accept$dataset == dataset_id & accept$role == role_full, , drop = FALSE]
    hand <- match(norm_file_key(files), norm_file_key(acc$file))
    for (i in which(!is.na(hand)))
      picks[[i]] <- pick(traces[[i]]$rt, traces[[i]]$it, acc$apex_rt[hand[i]] - APEX_PAD,
                         acc$apex_rt[hand[i]] + APEX_PAD, is_hilic)
    picks  <- add_msquality(chr, picks)
    keep   <- which(!vapply(picks, is.null, logical(1)))
    cols   <- paste0(role, "_", ROLE_COLS)
    if (!length(keep)) {
      m <- data.frame(file = character(0), stringsAsFactors = FALSE)
      for (cc in cols) m[[cc]] <- if (endsWith(cc, "detected") || endsWith(cc, "truncated"))
        logical(0) else if (endsWith(cc, "gate")) character(0) else numeric(0)
      return(m)
    }
    pk <- picks[keep]; fl <- files[keep]
    why <- vapply(seq_along(pk), function(i) {
      g <- if (is_dropped(dataset_id, role_full, fl[i], drop_all)) "manual drop"
           else gate_reason(pk[[i]], is_hilic)
      if (g == "accepted" && !is.na(hand[keep[i]])) "manual accept" else g
    }, character(1))
    det <- why %in% c("accepted", "manual accept")
    trunc <- vapply(pk, function(p) isTRUE(p$truncated), logical(1))
    m <- data.frame(
      file       = fl,
      auc        = ifelse(trunc, NA_real_, vapply(pk, function(p) p$auc, numeric(1))),
      auc0       = ifelse(trunc, NA_real_, vapply(pk, function(p) p$auc0, numeric(1))),
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
      truncated  = trunc,
      gate       = why,
      stringsAsFactors = FALSE)
    stats::setNames(m, c("file", cols))
  }

  p <- names(chr)
  metrics <- lapply(p, function(x)
    role_metrics(chr[[x]], x, role_of_prefix(x), cfg[[paste0("rt_", x, "_win")]]))
  for (i in seq_along(p))
    metrics[[i]][[paste0(p[i], "_fwhm_outlier")]] <- flag_fwhm(metrics[[i]], paste0(p[i], "_fwhm"))

  c(stats::setNames(peaks, paste0("peaks_", p)),
    stats::setNames(chr, paste0("chr_", p, "_narrow")),
    stats::setNames(metrics, paste0("metrics_", p)),
    list(rt_anchor = RT_ANCHOR))
}

# ---- 5. sample-table assembly --------------------------------------------

# to_bare_name() lives in R/metadata.R next to norm_file_key() (which calls it); several scripts
# source metadata.R alone, so it can't depend on ms1.R.

# Un-detect a peak of a target with a declared isobar (compound.csv `isobar`) that sits nearer its reference
# target's anchor than its own; only when both anchors are certain and >= CHAN_SEP apart. Returns `metrics`.
resolve_isobars <- function(metrics, conf, dataset_id,
                            rt_anchor = NULL,
                            spec = if (exists("compound_spec")) compound_spec else NULL) {
  # The role whose retention time each target's isobar shares; NA where none is declared.
  ref <- if (is.null(spec)) character(0) else stats::setNames(isobar_reference(spec), spec$role)
  # Each compound's centre: the median of its window's file anchors (MS2 times fired late on a tail
  # would pull it away); caches without stored anchors fall back to the raw MS2 times.
  centre <- function(p) {
    q <- rt_anchor[[p]]$anchors
    if (is.null(q))
      q <- as.numeric(conf$rt[conf$dataset_id == dataset_id & conf$role == role_of_prefix(p)])
    q <- q[is.finite(q)]
    if (length(q)) stats::median(q) else NA_real_
  }
  drop <- function(m, p) {
    col <- paste0(p, "_apex_rt"); det <- paste0(p, "_detected")
    r <- ref[role_of_prefix(p)]
    if (is.na(r) || !nzchar(r) || !all(c(col, det) %in% names(m))) return(m)
    q   <- target_prefix(r)
    own <- centre(p); oth <- centre(q)
    # Both anchors must be resolved: a split anchor's median sits between two peaks and would
    # reject peaks right where the compound elutes.
    resolved <- !isTRUE(rt_anchor[[p]]$uncertain) && !isTRUE(rt_anchor[[q]]$uncertain) &&
      is.finite(own) && is.finite(oth) && abs(own - oth) >= CHAN_SEP
    if (!resolved) return(m)
    bad <- m[[det]] %in% TRUE & is.finite(m[[col]]) &
      abs(m[[col]] - own) >= abs(m[[col]] - oth)
    m[[det]][bad] <- FALSE
    gcol <- paste0(p, "_gate")
    if (gcol %in% names(m)) m[[gcol]][bad] <- "isobar rule"
    m
  }
  stats::setNames(lapply(names(metrics), function(p) drop(metrics[[p]], p)), names(metrics))
}

# `metrics` holds each target's metrics, named by prefix, in compound.csv order (resolve_isobars()
# returns it so).
build_sample_table <- function(metrics, files_base, meta, dedup = TRUE) {
  # dedup = TRUE: one row per subject (avoids pseudoreplication); FALSE: one row per file. Merge on
  # norm_file_key, never the bare name: some cache names carry the accession twice.
  p <- names(metrics)
  for (x in p) {
    metrics[[x]]$.k   <- norm_file_key(to_bare_name(metrics[[x]]$file, files_base))
    metrics[[x]]$file <- NULL
  }
  # One merge per further target, in target order: with two targets, merge(par, met).
  sample_tbl <- Reduce(function(a, b) merge(a, b, by = ".k", all = TRUE), metrics)
  sample_tbl <- merge(data.frame(file = files_base, .k = norm_file_key(files_base),
                                 stringsAsFactors = FALSE),
                      sample_tbl, by = ".k", all.x = TRUE)
  sample_tbl$.k <- NULL

  # Join metadata on norm_file_key, not the raw name: a source may spell the extension differently
  # (.mzXML in a submitter TSV, .mzML in the deposit).
  sample_tbl$.k <- norm_file_key(sample_tbl$file)
  meta          <- meta[!duplicated(norm_file_key(meta$file)), , drop = FALSE]
  meta$.k       <- norm_file_key(meta$file)
  meta$file     <- NULL                     # keep the deposit's spelling
  sample_tbl    <- merge(sample_tbl, meta, by = ".k", all.x = TRUE)
  sample_tbl$.k <- NULL

  # Detection is decided once in the picker (<prefix>_detected); files with no pick merge in as NA -> FALSE.
  det   <- lapply(p, function(x) sample_tbl[[paste0(x, "_detected")]] %in% TRUE)
  roles <- role_of_prefix(p)

  if (length(p) == 2) {
    # Two targets: both, <role>_only (the metabolite's is metab_only) or neither, and the log2
    # ratio of the first target's area over the second's.
    only <- sub("^metabolite_only$", "metab_only", paste0(roles, "_only"))
    a <- det[[1]]; b <- det[[2]]
    sample_tbl$detection_tier <-
      ifelse( a &  b, "both",
      ifelse( a & !b, only[1],
      ifelse(!a &  b, only[2],
                      "neither")))

    sample_tbl$log2_ratio <- ifelse(
      sample_tbl$detection_tier == "both",
      log2(sample_tbl[[paste0(p[1], "_auc")]] / sample_tbl[[paste0(p[2], "_auc")]]), NA_real_)
  } else {
    # Any other number of targets: all, none, or the targets detected joined by "+"
    # ("parent+sulfone"). A ratio needs a pair, so there is no log2_ratio.
    tier <- character(nrow(sample_tbl))
    for (i in seq_along(p)) tier[det[[i]]] <- paste0(tier[det[[i]]], "+", roles[i])
    tier <- sub("^[+]", "", tier)
    tier[Reduce(`&`, det)]  <- "all"
    tier[!Reduce(`|`, det)] <- "none"
    sample_tbl$detection_tier <- tier
  }

  # Per-file table (subject_id retained) when dedup is off.
  if (!dedup) return(as.data.frame(sample_tbl))

  # One file per subject, highest combined AUC; the key falls back to sample_name, then file. Missing
  # key columns read as NA: a failed metadata fetch returns only `file`.
  col_or_na <- function(nm) if (nm %in% names(sample_tbl)) as.character(sample_tbl[[nm]])
                            else rep(NA_character_, nrow(sample_tbl))
  sample_tbl$.dedup_key <- dplyr::coalesce(
    col_or_na("subject_id"),
    col_or_na("sample_name"),
    as.character(sample_tbl$file))
  # The combined AUC, added target by target (a + b). Not rowSums(): it adds in long double,
  # which can round differently and flip a near-tie.
  sample_tbl$.auc_sum <- Reduce(`+`, lapply(p, function(x)
    dplyr::coalesce(sample_tbl[[paste0(x, "_auc")]], 0)))

  sample_tbl |>
    dplyr::group_by(.dedup_key) |>
    dplyr::slice_max(order_by = .auc_sum, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-.dedup_key, -.auc_sum) |>
    as.data.frame()
}

# ---- 6. cross-dataset wrapper --------------------------------------------

# The picker's behaviour, as a number: raise it whenever a change can alter a detection or a metric,
# and cached deposits are re-picked from their stored EICs (repick_cached).
PICKER_VERSION <- 5L

# The optional hand-authored picker inputs: apex overrides, drop-list and accept-list. A role not in
# compound.csv would silently match nothing, so it stops.
curation_inputs <- function(roles = compound_spec$role) {
  cur <- list(
    ov = if (file.exists("curation/rt_apex_overrides.csv"))
      read.csv("curation/rt_apex_overrides.csv", stringsAsFactors = FALSE) else
      data.frame(dataset_id = character(), role = character(), ms2_rt = numeric(), apex_rt = numeric()),
    drop_all = if (file.exists("curation/ms1_drop_list.csv"))
      read.csv("curation/ms1_drop_list.csv", stringsAsFactors = FALSE) else
      data.frame(dataset = character(), role = character(), file = character(), reason = character()),
    accept = if (file.exists("curation/ms1_accept_list.csv"))
      read.csv("curation/ms1_accept_list.csv", stringsAsFactors = FALSE) else
      data.frame(dataset = character(), role = character(), file = character(),
                 apex_rt = numeric(), reason = character()))
  bad <- setdiff(c(cur$ov$role, cur$drop_all$role, cur$accept$role), roles)
  if (length(bad))
    stop("curation/: role ", paste(bad, collapse = ", "), " is not a target in compound.csv (",
         paste(roles, collapse = ", "), ")")
  cur
}

# A deposit's own hand corrections, without their free-text reasons. Stored with the cache:
# adding or removing a correction changes the picks even though the EICs are unchanged.
curation_key <- function(dataset_id, cur = curation_inputs()) {
  d <- cur$drop_all[cur$drop_all$dataset == dataset_id, , drop = FALSE]
  o <- cur$ov[cur$ov$dataset_id == dataset_id, , drop = FALSE]
  a <- cur$accept[cur$accept$dataset == dataset_id, , drop = FALSE]
  sort(c(sprintf("drop %s %s", d$role, d$file),
         sprintf("apex %s %s %s", o$role, o$ms2_rt, o$apex_rt),
         sprintf("accept %s %s %s", a$role, a$file, a$apex_rt)))
}

# Re-run process_eics() on a cached deposit's stored EICs, as a fresh extraction would, without
# reading raw files; only the picking-dependent parts of the cache change.
repick_cached <- function(res, conf_final, files, ds_meta_tbl) {
  is_hilic <- row_is_hilic(dataset_meta_row(ds_meta_tbl, res$dataset_id))
  cur  <- curation_inputs()
  p    <- target_prefix(compound_spec$role)
  proc <- process_eics(per_target(res, "chr"), res$chr_all, res$dataset_id, res$cfg,
                       conf_final, cur$ov, cur$drop_all, is_hilic, basename(files))
  # Names already in the cache keep their place in it.
  for (k in c("rt_anchor", paste0("chr_", p, "_narrow"), paste0("peaks_", p),
              paste0("metrics_", p))) res[[k]] <- proc[[k]]
  res$sample_tbl <- build_sample_table(per_target(proc, "metrics"),
                                       basename(files), load_metadata(res$dataset_id))
  res$picker <- PICKER_VERSION; res$is_hilic <- is_hilic
  res$anchors <- anchor_key(conf_final, res$dataset_id)
  res$curation <- curation_key(res$dataset_id, cur)
  res
}

# The confirmed hits a deposit's peaks were anchored on, stored with the cache: a change makes the
# stored picks stale even when the EICs are not.
anchor_key <- function(conf_final, dataset_id) {
  h <- conf_final[conf_final$dataset_id == dataset_id, , drop = FALSE]
  sort(paste(h$role, basename(as.character(h$source_file)), round(as.numeric(h$rt), 2)))
}

# One metadata row per dataset, preferring a real chromatography where the harvest left several:
# the HILIC test needs length-1 fields.
dataset_meta_row <- function(ds_meta_tbl, dataset_id) {
  row <- ds_meta_tbl[ds_meta_tbl$dataset_id == dataset_id, ]
  if (!nrow(row))
    stop("No row for ", dataset_id, " in dataset_metadata_table")
  if (nrow(row) > 1) {
    real <- !is.na(row$chromatography) & !(row$chromatography %in% c("Unknown", ""))
    row  <- row[if (any(real)) which(real)[1] else 1, , drop = FALSE]
  }
  row
}
row_is_hilic <- function(row) {
  chrom <- row$chromatography; !is.na(chrom) && toupper(chrom) == "HILIC"
}

# Cross-dataset entry point. Resolves deposit_id/assay/mz/RT range from the artifacts; caller passes dataset_id + three artifacts.
run_ms1_for_dataset <- function(dataset_id,
                                conf_final,
                                bio_files_per_ds,
                                ds_meta_tbl,
                                ppm           = PPM_DEFAULT,
                                rt_pad_s      = 5,
                                cfg           = NULL,
                                verbose       = TRUE) {
  row <- dataset_meta_row(ds_meta_tbl, dataset_id)
  deposit_id <- row$deposit_id
  assay      <- row$assay

  files <- bio_files_per_ds[[dataset_id]]
  # The ion, m/z and RT range: 05 passes step 4's (extraction_window()); without them they are
  # derived here from the confirmations on the extracted files. The window is set in process_eics().
  if (is.null(cfg))
    cfg <- derive_extraction_config(dataset_id, conf_final,
                                    rt_pad_s = rt_pad_s, bio_files = files)
  # The m/z tolerance the traces are read at, kept with the cache: 05 re-extracts an assay
  # whose tolerance changed (05 passes the one step 4 set from the assay's instrument).
  cfg$ppm <- ppm
  meta  <- load_metadata(dataset_id)
  # Picker inputs: HILIC vs RP gates, manual apex overrides, and the reviewer drop-list.
  is_hilic <- row_is_hilic(row)
  cur <- curation_inputs(); ov <- cur$ov; drop_all <- cur$drop_all

  if (verbose) message("[", dataset_id, "] loading ", length(files),
                       " MS1 files …")
  # For MetaboLights, pass the assay so same-basename pos/neg modes are disambiguated at fetch time.
  ms1 <- load_spectra_ms1(deposit_id, files, assay = assay)

  roles <- compound_spec$role
  p     <- target_prefix(roles)

  if (verbose) message("[", dataset_id, "] extracting EICs …")
  eics <- extract_eics(ms1, cfg, ppm, roles)

  if (verbose) message("[", dataset_id, "] picking peaks + metrics …")
  proc <- process_eics(per_target(eics, "chr", p), eics$chr_all, dataset_id, cfg,
                       conf_final, ov, drop_all, is_hilic, basename(files))

  # One line per target, labels padded to one width ("Parent    ", "Metabolite").
  label <- format(paste0(toupper(substring(roles, 1, 1)), substring(roles, 2)))
  for (i in seq_along(roles)) {
    conf_files <- cfg[[paste0("confirmed_files_", roles[i])]]
    if (length(conf_files))
      check_confirmed(basename(dataOrigin(eics[[paste0("chr_", p[i])]])),
                      proc[[paste0("peaks_", p[i])]], conf_files, label[i])
  }

  sample_tbl <- build_sample_table(per_target(proc, "metrics", p),
                                   basename(files), meta)

  # The per-assay cache. Per target, in compound.csv order: chr_<prefix>, chr_<prefix>_narrow,
  # peaks_<prefix>, metrics_<prefix> (chr_par, chr_met, chr_par_narrow, ...).
  c(list(
    dataset_id     = dataset_id,
    deposit_id     = deposit_id,
    assay          = assay,
    cfg            = cfg,
    # Which picker, gate, anchors and hand corrections produced the metrics. 05
    # re-picks a cached deposit from its stored EICs when any no longer matches.
    picker         = PICKER_VERSION,
    is_hilic       = is_hilic,
    anchors        = anchor_key(conf_final, dataset_id),
    curation       = curation_key(dataset_id, cur),
    # Per role: the window used, how many files gave an anchor, and whether the anchors split between
    # equally supported RTs (uncertain = TRUE: an unresolved anchor, not a confident wide peak).
    rt_anchor      = proc$rt_anchor,
    chr_all        = eics$chr_all),
    eics[paste0("chr_", p)],
    proc[c(paste0("chr_", p, "_narrow"), paste0("peaks_", p), paste0("metrics_", p))],
    list(sample_tbl = sample_tbl))
}

# One per-target field of an extraction result, as a list named by prefix in target order:
# per_target(r, "metrics") is list(par = r$metrics_par, met = r$metrics_met).
per_target <- function(x, field, prefixes = target_prefix(compound_spec$role)) {
  stats::setNames(lapply(prefixes, function(p) x[[paste0(field, "_", p)]]), prefixes)
}

# ===== peak detection + metrics =====

# Count peaks per EIC via a local-maxima detector. hws: half-window in scans (max within ±hws
# neighbours). pct_max: min height as a fraction of the trace max (filters sub-noise wiggles).
count_peaks <- function(chr, hws = 3L, pct_max = 0.1) {
  ints <- intensity(chr)
  vapply(ints, function(x) {
    x_s <- x[!is.na(x)]
    if (!length(x_s)) return(0L)   # an empty trace has no peaks; max() on it only warns
    mx <- max(x_s)
    if (!is.finite(mx)) return(0L)
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
  # norm_file_key(), not a private copy of the rule, so the two cannot drift apart.
  idx <- which(norm_file_key(files) %in% norm_file_key(conf_files))
  cat(role, "- confirmed files loaded:", length(idx),
      "| with >=1 peak:", sum(n_peaks[idx] > 0), "\n")

  missed <- files[idx[n_peaks[idx] == 0]]
  if (length(missed))
    cat("  WARNING missed:", paste(missed, collapse = ", "), "\n")
}

