## download_confirmed_massive.R
## Download mzML/mzXML files for MassIVE datasets that passed SIRIUS
## confirmation (strict or standard), filtering out QC/blank/pool samples
## using both SampleType and filename columns from dataset_metadata/.
## Run on cluster node:
##   Rscript download_confirmed_massive.R

suppressPackageStartupMessages(library(MsBackendMassIVE))

# ── patterns that identify non-biological files ──────────────────────────────
EXCLUDE_SAMPLETYPE_RE <- "(?i)^(blank|qc|pool|reference|standard|control)"
EXCLUDE_FILENAME_RE   <- "(?i)(^|[_\\-/])(blank|QC|pool|reference|standard|control|mix)[_\\-\\.]"

# ── helper: load one metadata file (TSV or CSV) ───────────────────────────────
load_metadata <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "tsv") {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# ── helper: filter to biological samples only ─────────────────────────────────
filter_biological <- function(meta, ds) {
  n_start <- nrow(meta)

  # 1. SampleType filter (when column is present)
  if ("SampleType" %in% names(meta)) {
    keep_st <- !grepl(EXCLUDE_SAMPLETYPE_RE, meta$SampleType, perl = TRUE)
    n_excl  <- sum(!keep_st)
    if (n_excl > 0)
      cat("  [", ds, "] SampleType filter: dropped", n_excl, "rows\n")
    meta <- meta[keep_st, ]
  } else {
    cat("  [", ds, "] no SampleType column — skipping SampleType filter\n")
  }

  # 2. filename filter (always applied)
  if ("filename" %in% names(meta)) {
    keep_fn <- !grepl(EXCLUDE_FILENAME_RE, basename(meta$filename), perl = TRUE)
    n_excl  <- sum(!keep_fn)
    if (n_excl > 0)
      cat("  [", ds, "] filename filter:   dropped", n_excl, "rows\n")
    meta <- meta[keep_fn, ]
  } else {
    warning("No 'filename' column found in metadata for ", ds)
  }

  cat("  [", ds, "] kept", nrow(meta), "of", n_start, "files\n")
  meta
}

# ── main ──────────────────────────────────────────────────────────────────────
confirmed <- read.csv("hits-confirmed-massive.csv", check.names = FALSE)
datasets  <- unique(confirmed$dataset_id)
cat("Datasets to download:", length(datasets), "\n")
cat(paste(" ", datasets, collapse = "\n"), "\n\n")

meta_dir <- "dataset_metadata"

for (ds in datasets) {
  cat("[", ds, "] resolving files to download ...\n")

  # locate metadata file for this dataset (TSV preferred over CSV)
  meta_files <- list.files(meta_dir,
                           pattern  = paste0("^", ds, "_metadata\\.(tsv|csv)$"),
                           full.names = TRUE)

  wanted <- NULL

  if (length(meta_files) == 0) {
    cat("  [", ds, "] WARNING: no metadata file found — downloading all files\n")
  } else {
    if (length(meta_files) > 1) {
      # prefer TSV
      tsv <- meta_files[grepl("\\.tsv$", meta_files)]
      meta_files <- if (length(tsv)) tsv[1] else meta_files[1]
    }
    cat("  [", ds, "] using metadata:", basename(meta_files), "\n")
    meta    <- load_metadata(meta_files)
    meta    <- filter_biological(meta, ds)
    wanted  <- if ("filename" %in% names(meta)) unique(meta$filename) else NULL
  }

  # download
  cat("[", ds, "] starting download ...\n")
  tryCatch({
    if (is.null(wanted)) {
      massive_sync_data_files(massiveId = ds)
    } else {
      cat("  [", ds, "] files to download:", length(wanted), "\n")
      massive_sync_data_files(massiveId = ds, fileName = wanted)
    }
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

  cat("[", ds, "] done\n\n")
}

cat("All downloads complete.\n")
