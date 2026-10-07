# Rscript R/peak_review_panels.R <application> [dataset ...]: re-picks the peaks with R/ms1.R's
# picker from cached wide EICs; writes panels + index.csv to application/<name>/peaks_review_v2/.
ROOT <- dirname(dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1]))))
args <- commandArgs(trailingOnly = TRUE)
APP  <- args[1]
if (is.na(APP) || !dir.exists(file.path(ROOT, "application", APP)))
  stop("usage: Rscript R/peak_review_panels.R <application> [dataset ...]")
setwd(file.path(ROOT, "application", APP))   # first: setup.R creates artifacts/ in cwd
suppressPackageStartupMessages({
  source(file.path(ROOT, "R/setup.R")); source(file.path(ROOT, "R/metadata.R"))
  source(file.path(ROOT, "R/ms1.R")); library(Chromatograms)
})
# Single-threaded: the cost is disk I/O, so parallel backends do not help.
BiocParallel::register(BiocParallel::SerialParam())

# Picker constants + functions (walk_bound, pick, accept_peak, acceptance_window, snap_one, is_dropped)
# all come from R/ms1.R (sourced above). Only the panel-runner extras live here.
PLOT_PAD <- 4                      # RT margin (s) around the peak in the panels
OUT      <- "peaks_review_v2"
dir.create(OUT, showWarnings = FALSE)

bfd <- readRDS("artifacts/bio_files_per_ds.rds")
dm  <- read.csv("artifacts/dataset_metadata_table.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Confirmed-hit anchors + manual overrides + drop-list, fed to R/ms1.R's acceptance_window() / is_dropped().
hits_all <- read.csv("artifacts/hits-confirmed-all.csv", stringsAsFactors = FALSE)
ov_all <- if (file.exists("curation/rt_apex_overrides.csv"))
  read.csv("curation/rt_apex_overrides.csv", stringsAsFactors = FALSE) else
  data.frame(dataset_id = character(), role = character(), ms2_rt = numeric(), apex_rt = numeric())
drop_all <- if (file.exists("curation/ms1_drop_list.csv"))
  read.csv("curation/ms1_drop_list.csv", stringsAsFactors = FALSE) else
  data.frame(dataset = character(), role = character(), file = character(), reason = character())
accept_all <- if (file.exists("curation/ms1_accept_list.csv"))
  read.csv("curation/ms1_accept_list.csv", stringsAsFactors = FALSE) else
  data.frame(dataset = character(), role = character(), file = character(), apex_rt = numeric())

wide_eic <- function(chr_all, mz, ppm, rt_lo, rt_hi) {
  orig <- unique(dataOrigin(chr_all))
  pt <- data.frame(msLevel = 1L, dataOrigin = orig, rtMin = rt_lo, rtMax = rt_hi,
                   mzMin = mz * (1 - ppm * 1e-6), mzMax = mz * (1 + ppm * 1e-6),
                   stringsAsFactors = FALSE)
  chromExtract(chr_all, pt, by = c("msLevel", "dataOrigin"))
}

panel <- function(path, p, win, tag, fname, sub, det, raw = NULL) {
  png(path, width = 900, height = 560, res = 108); on.exit(dev.off())
  ctr <- mean(win)                                   # expected RT: always keep it in frame
  xlo <- min(p$lb, ctr) - PLOT_PAD; xhi <- max(p$rb, ctr) + PLOT_PAD   # boundaries already span broad peaks
  keep <- p$rt >= xlo & p$rt <= xhi
  if (!any(keep)) { plot.new(); title(tag); return() }
  plot(NA, xlim = c(xlo, xhi), ylim = c(0, max(p$it[keep]) * 1.08),
       xlab = "RT (s)", ylab = "Intensity", main = tag, cex.main = 0.82,
       col.main = if (det) "#1B7837" else "#B2182B")
  rect(win[1], 0, win[2], max(p$it[keep]) * 1.08, col = adjustcolor("#4393C3", 0.10), border = NA)
  abline(v = mean(win), col = "#2166AC", lty = 2, lwd = 1.1)   # expected RT (window centre)
  pk <- p$idx
  bl <- if (p$rb > p$lb) p$lb_int + (p$rb_int - p$lb_int) * (p$rt[pk] - p$lb) / (p$rb - p$lb)
        else rep(min(p$lb_int, p$rb_int), length(pk))       # local baseline line
  polygon(c(p$rt[pk], rev(p$rt[pk])), c(p$it[pk], rev(bl)),  # shade peak ABOVE the baseline
          col = adjustcolor("#D6604D", 0.22), border = NA)
  lines(p$rt[pk], bl, col = "#D6604D", lty = 3, lwd = 0.9)   # the subtracted baseline
  lines(p$rt[keep], p$it[keep], col = "grey45", lwd = 1.3)
  points(p$rt[keep], p$it[keep], pch = 16, cex = 0.7, col = "grey25")
  # Empty scans (nothing inside the mass window) drawn on the axis: the picker drops them, and
  # without them a peak cut off by its own m/z window looks like one sitting on a baseline.
  if (!is.null(raw)) {
    z <- raw$rt[!is.finite(raw$it) | raw$it <= 0]; z <- z[z >= xlo & z <= xhi]
    points(z, rep(0, length(z)), pch = 1, cex = 0.8, col = "grey40")
  }
  points(p$rt[pk], p$it[pk], pch = 16, cex = 0.9, col = "#D6604D")
  abline(v = p$apex_rt, col = "#1B7837", lty = 3)
  mtext(fname, side = 3, line = 0.95, cex = 0.66, col = "grey40")
  mtext(sub,   side = 3, line = 0.15, cex = 0.70, col = "grey30")
}

sel <- args[-1]
if (!length(sel)) sel <- names(sort(vapply(bfd, length, integer(1))))

dir.create(file.path(OUT, "_cache"), showWarnings = FALSE, recursive = TRUE)

# Caches the raw (rt, it) traces of the wide EICs, one entry per target prefix; rebuilt when
# the assay's m/z tolerance changes or a target is missing.
build_cache <- function(ds, r) {
  ppm <- r$cfg$ppm
  P   <- names(r$rt_anchor)
  cf <- file.path(OUT, "_cache", paste0(ds, ".rds"))
  if (file.exists(cf)) {
    cache <- readRDS(cf)
    if (isTRUE(cache$ppm == ppm) && all(P %in% names(cache))) return(cache)
  }
  t0 <- proc.time()[3]
  cache <- list(ppm = ppm)
  for (role in P) {
    win <- r$cfg[[paste0("rt_", role, "_win")]]
    mz  <- r$cfg[[paste0("mz_", role)]]
    chr <- Chromatograms::filterEmptyChromatograms(
             wide_eic(r$chr_all, mz, ppm, win[1] - WIDE_EXTRA, win[2] + WIDE_EXTRA))
    pd  <- if (length(chr)) Chromatograms::peaksData(chr) else list()
    cache[[role]] <- list(
      win = win, files = basename(dataOrigin(chr)),
      traces = lapply(pd, function(d) list(rt = d[, "rtime"], it = d[, "intensity"])))
  }
  secs <- round(proc.time()[3] - t0, 1)
  nf <- length(unique(unlist(lapply(P, function(role) cache[[role]]$files))))
  tl <- file.path(OUT, "_extract_times.csv")
  # Stamp the run date and the Chromatograms version: extraction cost tracks that
  # package, so a timing is only comparable to another run on the same version.
  if (!file.exists(tl))
    cat("dataset,n_files,extract_seconds,run_date,chromatograms\n", file = tl)
  cat(sprintf("%s,%d,%.1f,%s,%s\n", ds, nf, secs, format(Sys.Date()),
              as.character(utils::packageVersion("Chromatograms"))),
      file = tl, append = TRUE)
  message(sprintf("[%s] EIC extraction: %.1f s for %d files", ds, secs, nf))
  saveRDS(cache, cf); cache
}

for (ds in sel) {
  f <- file.path("artifacts/per_dataset", paste0(ds, ".rds")); if (!file.exists(f)) next
  r <- readRDS(f)
  if (is.null(r$cfg$ppm)) {
    message("[", ds, "] skipped: its step 5 cache records no m/z tolerance; render 05 first")
    next
  }
  chrom <- dm$chromatography[match(ds, dm$dataset_id)]
  is_hilic <- !is.na(chrom) && toupper(chrom) == "HILIC"
  fb <- basename(as.character(bfd[[ds]]))
  bare <- function(x) norm_file_key(to_bare_name(basename(as.character(x)), fb))
  dd <- file.path(OUT, ds); dir.create(dd, showWarnings = FALSE)
  for (s in c("detected", "filtered")) unlink(file.path(dd, s), recursive = TRUE)
  for (s in c("detected", "filtered")) dir.create(file.path(dd, s), showWarnings = FALSE)
  cache <- build_cache(ds, r)
  P <- names(r$rt_anchor)   # the targets' prefixes (par, met), in target order
  rows <- list(); n_any_auc <- setNames(numeric(length(P)), P); n_accepted <- n_any_auc
  for (role in P) {
    role_full <- role_of_prefix(role)
    win <- acceptance_window(ds, role_full, cache[[role]], bare, cache[[role]]$win, hits_all, ov_all, is_hilic)
    bn <- cache[[role]]$files; tr <- cache[[role]]$traces
    m <- as.data.frame(r[[paste0("metrics_", role)]])
    mk <- basename(m$file)
    for (i in seq_along(tr)) {
      # A peak accepted by hand is picked around the reviewer's apex, as in step 5.
      a <- accept_all[accept_all$dataset == ds & accept_all$role == role_full &
                      norm_file_key(accept_all$file) == norm_file_key(bn[i]), , drop = FALSE]
      p <- if (nrow(a)) pick(tr[[i]]$rt, tr[[i]]$it, a$apex_rt[1] - APEX_PAD,
                             a$apex_rt[1] + APEX_PAD, is_hilic)
           else pick(tr[[i]]$rt, tr[[i]]$it, win[1], win[2], is_hilic)
      id <- sprintf("%s_%s_%03d", ds, role, i)
      dropped <- is_dropped(ds, role_full, bn[i], drop_all)   # manual reviewer removal
      why <- if (dropped) "manual drop" else gate_reason(p, is_hilic)
      if (nrow(a) && why == "accepted") why <- "manual accept"
      det <- why %in% c("accepted", "manual accept")
      n_accepted[role] <- n_accepted[role] + det
      j <- match(bn[i], mk)
      any_auc <- if (!is.na(j)) m[[paste0(role, "_auc")]][j] else NA
      if (!is.na(any_auc)) n_any_auc[role] <- n_any_auc[role] + 1  # any integrated peak, gate aside
      if (is.null(p)) {   # nothing to draw, but still a rejection to count
        rows[[length(rows) + 1L]] <- data.frame(id, role, file = bn[i], detected = FALSE,
          reason = why, truncated = NA, pts = NA, fwhm = NA, snr = NA, pkwidth = NA, auc = NA, apex_rt = NA,
          png = NA, verdict = "", notes = "", stringsAsFactors = FALSE)
        next
      }
      sd <- if (det) "detected" else "filtered"
      tag <- sprintf("%s | %s | %s", id, role_full,
                     if (det) "DETECTED" else if (dropped) "MANUAL DROP" else paste("filtered:", why))
      sub <- sprintf("FWHM %s s | SNR %s | pkwidth %.1f s | pts %d | AUC %s%s",
                     if (is.na(p$fwhm)) "NA" else round(p$fwhm, 1),
                     if (is.finite(p$snr)) round(p$snr, 1) else "NA",
                     p$pkwidth, p$npts, signif(p$auc, 3),
                     if (isTRUE(p$truncated)) " | truncated: area withheld" else "")
      panel(file.path(dd, sd, paste0(id, ".png")), p, win, tag, bn[i], sub, det, raw = tr[[i]])
      rows[[length(rows) + 1L]] <- data.frame(id, role, file = bn[i], detected = det,
        reason = why, truncated = isTRUE(p$truncated), pts = p$npts, fwhm = p$fwhm, snr = p$snr,
        pkwidth = p$pkwidth, auc = p$auc,
        apex_rt = p$apex_rt, png = file.path(ds, sd, paste0(id, ".png")),
        verdict = "", notes = "", stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) write.csv(do.call(rbind, rows), file.path(dd, "index.csv"), row.names = FALSE)
  message(ds, ": ", paste(sprintf("%s any-AUC=%d accepted=%d", role_of_prefix(P),
                                  n_any_auc, n_accepted), collapse = " | "))
}
message("done -> ", OUT, "/")
