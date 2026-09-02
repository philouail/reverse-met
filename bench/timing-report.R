# Print the pipeline timing table: notebook renders (artifacts/pipeline_timing.csv,
# from bench/time-pipeline.sh) and MS1 extraction (peaks_review_v2/_extract_times.csv,
# from R/peak_review_panels.R).
#   Rscript bench/timing-report.R
setwd(dirname(dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1])))))

fmt <- function(s) if (is.na(s)) "-" else if (s < 90) sprintf("%.0f s", s) else
  if (s < 5400) sprintf("%.1f min", s / 60) else sprintf("%.2f h", s / 3600)

cat("\n===== ANALYSIS NOTEBOOKS =====\n")
tp <- "artifacts/pipeline_timing.csv"
if (file.exists(tp)) {
  cat(paste(grep("^#", readLines(tp), value = TRUE), collapse = "\n"), "\n\n")
  d <- read.csv(tp, comment.char = "#", stringsAsFactors = FALSE)
  d$time <- vapply(d$seconds, fmt, character(1))
  d$ok   <- ifelse(d$exit == 0, "ok", paste0("FAILED (", d$exit, ")"))
  print(d[, c("stage", "notebook", "time", "ok")], row.names = FALSE)
  cat("\ntotal:", fmt(sum(d$seconds)), "\n")
} else cat("(no pipeline_timing.csv)\n")

cat("\n===== MS1 EXTRACTION, per deposit =====\n")
ep <- "peaks_review_v2/_extract_times.csv"
if (file.exists(ep)) {
  e <- read.csv(ep, stringsAsFactors = FALSE)
  e <- e[!duplicated(e$dataset, fromLast = TRUE), ]
  e <- e[order(-e$n_files), ]
  e$time <- vapply(e$extract_seconds, fmt, character(1))
  e$s_per_file <- round(e$extract_seconds / e$n_files, 2)
  cols <- c("dataset", "n_files", "time", "s_per_file",
            intersect(c("run_date", "chromatograms"), names(e)))
  print(e[, cols], row.names = FALSE)
  cat(sprintf("\n%d deposits, %s files, %s, %.2f s per file\n",
              nrow(e), format(sum(e$n_files), big.mark = ","),
              fmt(sum(e$extract_seconds)),
              sum(e$extract_seconds) / sum(e$n_files)))
} else cat("(no _extract_times.csv)\n")
cat("\n")
