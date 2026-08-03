# Peak-review panel generator for the manual audit. Uses the CANONICAL picker from R/ms1.R
# (constants + walk_bound + pick + detected_new + new_win + snap_one + is_dropped, sourced below);
# this file adds only the panel-drawing + wide-EIC cache runner. Does NOT touch the pipeline.
# Re-picks the confirmation peaks from the cached wide EIC (chr_all): builds the MS1-apex window
# (new_win), integrates each in-window local max valley-to-valley (pick), gates it (detected_new),
# applies the drop-list; writes panels to peaks_review_v2/<ds>/{detected,filtered} + index.csv.
#   Rscript R/peak_review_panels.R MSV000094097 [more datasets...]
suppressPackageStartupMessages({
  source("R/_setup.R"); source("R/metadata.R"); source("R/ms1.R"); library(Chromatograms)
})
# Single-threaded extraction: the parallel chromExtract spawns a large worker swarm
# that thrashes the machine and is SLOWER here; serial is lighter.
BiocParallel::register(BiocParallel::SerialParam())

# Picker constants + functions (walk_bound, pick, detected_new, new_win, snap_one, is_dropped)
# all come from R/ms1.R (sourced above). Only the panel-runner extras live here.
PLOT_PAD <- 4                      # RT margin each side of the peak boundaries in review panels
OUT      <- "peaks_review_v2"
dir.create(OUT, showWarnings = FALSE)

bfd <- readRDS("artifacts/bio_files_per_ds.rds")
dm  <- read.csv("artifacts/dataset_metadata_table.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Confirmed-hit anchors + manual overrides + drop-list, fed to R/ms1.R's new_win() / is_dropped().
hits_all <- read.csv("artifacts/hits-confirmed-all.csv", stringsAsFactors = FALSE)
ov_all <- if (file.exists("curation/rt_apex_overrides.csv"))
  read.csv("curation/rt_apex_overrides.csv", stringsAsFactors = FALSE) else
  data.frame(dataset_id = character(), role = character(), ms2_rt = numeric(), apex_rt = numeric())
drop_all <- if (file.exists("curation/ms1_drop_list.csv"))
  read.csv("curation/ms1_drop_list.csv", stringsAsFactors = FALSE) else
  data.frame(dataset = character(), role = character(), file = character(), reason = character())

wide_eic <- function(chr_all, mz, ppm, rt_lo, rt_hi) {
  orig <- unique(dataOrigin(chr_all))
  pt <- data.frame(msLevel = 1L, dataOrigin = orig, rtMin = rt_lo, rtMax = rt_hi,
                   mzMin = mz * (1 - ppm * 1e-6), mzMax = mz * (1 + ppm * 1e-6),
                   stringsAsFactors = FALSE)
  chromExtract(chr_all, pt, by = c("msLevel", "dataOrigin"))
}

panel <- function(path, p, win, tag, fname, sub, det) {
  png(path, width = 900, height = 560, res = 108); on.exit(dev.off())
  ctr <- mean(win)                                   # expected RT — always keep it in frame
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
  points(p$rt[pk], p$it[pk], pch = 16, cex = 0.9, col = "#D6604D")
  abline(v = p$apex_rt, col = "#1B7837", lty = 3)
  mtext(fname, side = 3, line = 0.95, cex = 0.66, col = "grey40")
  mtext(sub,   side = 3, line = 0.15, cex = 0.70, col = "grey30")
}

sel <- commandArgs(trailingOnly = TRUE)
if (!length(sel)) sel <- names(sort(vapply(bfd, length, integer(1))))

dir.create(file.path(OUT, "_cache"), showWarnings = FALSE, recursive = TRUE)

# Slow one-time step: pull the wide EICs off the Spectra-backed chr_all and store
# the raw (rt,it) traces so panel re-runs are instant.
build_cache <- function(ds, r, ppm = 20) {
  short <- sub("_(LC-MS|AN[0-9]).*$", "", ds)
  cf <- file.path(OUT, "_cache", paste0(short, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  t0 <- proc.time()[3]
  cache <- list()
  for (role in c("par", "met")) {
    win <- if (role == "par") r$cfg$rt_par_win else r$cfg$rt_met_win
    mz  <- if (role == "par") r$cfg$mz_par     else r$cfg$mz_met
    chr <- Chromatograms::filterEmptyChromatograms(
             wide_eic(r$chr_all, mz, ppm, win[1] - WIDE_EXTRA, win[2] + WIDE_EXTRA))
    pd  <- if (length(chr)) Chromatograms::peaksData(chr) else list()
    cache[[role]] <- list(
      win = win, files = basename(dataOrigin(chr)),
      traces = lapply(pd, function(d) list(rt = d[, "rtime"], it = d[, "intensity"])))
  }
  secs <- round(proc.time()[3] - t0, 1)
  nf <- length(unique(c(cache$par$files, cache$met$files)))
  tl <- file.path(OUT, "_extract_times.csv")
  if (!file.exists(tl)) cat("dataset,n_files,extract_seconds\n", file = tl)  # header on first write
  cat(sprintf("%s,%d,%.1f\n", short, nf, secs), file = tl, append = TRUE)
  message(sprintf("[%s] EIC extraction: %.1f s for %d files", short, secs, nf))
  saveRDS(cache, cf); cache
}

for (ds in sel) {
  f <- file.path("artifacts/per_dataset", paste0(ds, ".rds")); if (!file.exists(f)) next
  r <- readRDS(f); short <- sub("_(LC-MS|AN[0-9]).*$", "", ds)
  chrom <- dm$chromatography[match(ds, dm$dataset_id)]
  is_hilic <- !is.na(chrom) && toupper(chrom) == "HILIC"
  fb <- basename(as.character(bfd[[ds]]))
  bare <- function(x) norm_file_key(to_bare_name(basename(as.character(x)), fb))
  dd <- file.path(OUT, short); dir.create(dd, showWarnings = FALSE)
  for (s in c("detected", "filtered")) unlink(file.path(dd, s), recursive = TRUE)
  for (s in c("detected", "filtered")) dir.create(file.path(dd, s), showWarnings = FALSE)
  cache <- build_cache(ds, r)
  rows <- list(); nd_old <- c(par = 0, met = 0); nd_new <- c(par = 0, met = 0)
  for (role in c("par", "met")) {
    role_full <- if (role == "par") "parent" else "metabolite"
    win <- new_win(ds, role_full, cache[[role]], bare, cache[[role]]$win, hits_all, ov_all)
    bn <- cache[[role]]$files; tr <- cache[[role]]$traces
    m <- as.data.frame(if (role == "par") r$metrics_par else r$metrics_met)
    mk <- basename(m$file)
    for (i in seq_along(tr)) {
      p <- pick(tr[[i]]$rt, tr[[i]]$it, win[1], win[2], is_hilic)
      id <- sprintf("%s_%s_%03d", short, role, i)
      dropped <- is_dropped(short, role_full, bn[i], drop_all)   # manual reviewer removal
      det <- detected_new(p, is_hilic) && !dropped
      nd_new[role] <- nd_new[role] + det
      j <- match(bn[i], mk)
      old_auc <- if (!is.na(j)) m[[paste0(role, "_auc")]][j] else NA
      if (!is.na(old_auc)) nd_old[role] <- nd_old[role] + 1  # old counted any AUC
      if (is.null(p)) next
      sd <- if (det) "detected" else "filtered"
      tag <- sprintf("%s | %s | %s", id, if (role == "par") "parent" else "metabolite",
                     if (det) "DETECTED (new)" else if (dropped) "MANUAL DROP" else "filtered (new)")
      sub <- sprintf("FWHM %s s | SNR %s | pkwidth %.1f s | pts %d | AUC %s",
                     if (is.na(p$fwhm)) "NA" else round(p$fwhm, 1),
                     if (is.finite(p$snr)) round(p$snr, 1) else "NA",
                     p$pkwidth, p$npts, signif(p$auc, 3))
      panel(file.path(dd, sd, paste0(id, ".png")), p, win, tag, bn[i], sub, det)
      rows[[length(rows) + 1L]] <- data.frame(id, role, file = bn[i], detected = det,
        pts = p$npts, fwhm = p$fwhm, snr = p$snr, pkwidth = p$pkwidth, auc = p$auc,
        apex_rt = p$apex_rt, png = file.path(short, sd, paste0(id, ".png")),
        verdict = "", notes = "", stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) write.csv(do.call(rbind, rows), file.path(dd, "index.csv"), row.names = FALSE)
  message(sprintf("%s: parent old(any-AUC)=%d new=%d | metab old=%d new=%d",
                  short, nd_old["par"], nd_new["par"], nd_old["met"], nd_new["met"]))
}
message("done -> ", OUT, "/")
