# R/metadata.R — per-dataset biological metadata, table-driven.
# harvest_metadata_coverage(ds): auto-match canonical -> native column names.
# load_metadata(dataset_id): apply the table's col_*/const_* spec to the native table.

CANONICAL_COLS <- c("file", "sample_name", "subject_id", "sex", "age",
                    "group", "health_status", "body_site", "country",
                    "treatment", "timepoint")

# MIN_SUBMITTER_FIELDS: min canonical DRAFT_CANDIDATES fields a submitter-only file
# must populate to be admitted (see msv_native_sample_table).
# MSV_FILENAME_COLS: native raw-file-name columns; naming varies across submitter TSVs.
MIN_SUBMITTER_FIELDS <- 2L

MSV_FILENAME_COLS <- c("filename", "Filename", "FileName",
                       "Metabolomics_FileName")

# Combined dataset_id key: one per (deposit, assay); bare accession when assay is NA.
make_dataset_id <- function(deposit_id, assay) {
  mapply(function(d, a) {
    if (is.na(a) || !nzchar(a)) return(d)
    bare <- sub("\\.(txt|csv|tsv)$", "", as.character(a))
    bare <- sub(paste0("^a_", d, "_?"), "", bare)
    bare <- sub(paste0("^",   d, "_?"), "", bare)
    if (!nzchar(bare)) d else paste0(d, "_", bare)
  }, deposit_id, assay, USE.NAMES = FALSE)
}

# Candidate native columns per canonical field; matched case-insensitively full-string,
# first hit wins. Covers ISA-Tab (MetaboLights), mwTab (Workbench), MassIVE "ATTRIBUTE_*",
# and Pan-ReDU bare names.
DRAFT_CANDIDATES <- list(
  subject_id    = c("Source Name", "Subject", "Subject ID", "SubjectID",
                    "ATTRIBUTE_Subject",
                    "ATTRIBUTE_SubjectIdentifierAsRecorded",
                    "SubjectIdentifierAsRecorded", "UniqueSubjectID"),
  sample_name   = c("Sample Name", "SampleID", "Sample ID",
                    "ATTRIBUTE_SampleName", "qiita_sample_name", "filename"),
  group         = c("Comment\\[Patient status\\]",
                    "Factor Value\\[.*\\]",
                    "Factors: Phenotype", "Factors: Disease.*",
                    "Factors: Health.*", "Factors: Status.*",
                    "Factors: Group.*", "Factors: Treatment.*",
                    "Comment\\[Group\\]",
                    "ATTRIBUTE_Disease", "ATTRIBUTE_HealthStatus",
                    "ATTRIBUTE_Group",
                    "HealthStatus", "DOIDCommonName"),
  health_status = c("Comment\\[Patient status\\]",
                    "Factor Value\\[Disease\\]",
                    "Factor Value\\[Health.*\\]",
                    "Factors: Phenotype", "Factors: Disease.*",
                    "Factors: Health.*", "Factors: Status.*",
                    "ATTRIBUTE_HealthStatus", "ATTRIBUTE_Disease",
                    "HealthStatus"),
  # Keep both BodySite spellings: 082493's underscored ATTRIBUTE_Body_Site else NA-s
  # its 341 blood_plasma files, the matrix the CYP2C19 ratio is defined in.
  body_site     = c("Characteristics\\[Organism part\\]",
                    "Characteristics\\[Sample type\\]",
                    "Characteristics\\[Body site\\]",
                    "ATTRIBUTE_BodySite", "ATTRIBUTE_Body_Site",
                    "ATTRIBUTE_Sampletype", "ATTRIBUTE_Sample_Type",
                    "UBERONBodyPartName"),
  # Sampling time: needed to reconstruct a PK curve and integrate per-subject AUC.
  timepoint     = c("ATTRIBUTE_Time_Point_Mins", "ATTRIBUTE_Time_Point",
                    "ATTRIBUTE_Timepoint", "ATTRIBUTE_Study_Day",
                    "Factor Value\\[Time.*\\]", "Factors: Time.*",
                    "Characteristics\\[Time.*\\]", "collection_timestamp"),
  sex           = c("Characteristics\\[Sex\\]",
                    "Factor Value\\[Sex\\]",
                    "Comment\\[Sex\\]",
                    "Comment\\[Gender\\]",
                    "Factors: Gender", "Factors: Sex",
                    "ATTRIBUTE_Sex",
                    "BiologicalSex"),
  age           = c("Characteristics\\[Age\\]",
                    "Factor Value\\[Age\\]",
                    "Comment\\[Age\\]",
                    "Additional sample data: Age",
                    "ATTRIBUTE_Age",
                    "AgeInYears"),
  country       = c("Comment\\[Country\\]",
                    "Characteristics\\[Geographic location\\]",
                    "Characteristics\\[Country\\]",
                    "ATTRIBUTE_Country",
                    "Country"),
  treatment     = c("Factor Value\\[Treatment\\]",
                    "Comment\\[Treatment\\]",
                    "Characteristics\\[Treatment\\]",
                    "ATTRIBUTE_Treatment")
)

match_candidate <- function(native_cols, patterns) {
  for (pat in patterns) {
    hits <- grep(paste0("^", pat, "$"), native_cols,
                 value = TRUE, ignore.case = TRUE)
    if (length(hits)) return(hits[1])
  }
  NA_character_
}

# ---- harvest --------------------------------------------------------------

native_sample_table <- function(deposit_id) {
  if (startsWith(deposit_id, "MTBLS")) {
    # Cached ISA-Tab study table first, same as the assay table in
    # expand_files_mtbls(): a render must not depend on EBI being reachable.
    cache <- file.path(META_DIR, paste0("s_", deposit_id, ".txt"))
    if (file.exists(cache)) read.delim(cache, check.names = FALSE,
                                       stringsAsFactors = FALSE) else
      as.data.frame(MsBackendMetaboLights::mtbls_sample_data(deposit_id),
                    check.names = FALSE)
  } else if (startsWith(deposit_id, "ST")) {
    MsBackendMetabolomicsWorkbench::mwb_metadata(deposit_id)$sample_annotation
  } else if (startsWith(deposit_id, "MSV")) {
    msv_native_sample_table(deposit_id)
  } else NULL
}

# MassIVE: full outer join of Pan-ReDU and submitter TSV on filename stem, then drop
# submitter-only rows absent from massive_list_files(). Pan-ReDU under-covers deposits
# (misses 082493's 341-file blood collection + its CYP2C19 columns); the inventory check
# keeps phantom protection while recovering annotated files. NULL if both sources missing.
msv_native_sample_table <- function(deposit_id) {
  pr <- fetch_panredu_metadata(deposit_id)

  fs <- list.files("dataset_metadata",
                   pattern = paste0("^", deposit_id, "_metadata\\.(csv|tsv)$"),
                   full.names = TRUE)
  local <- if (length(fs)) {
    sep <- if (endsWith(fs[1], ".tsv")) "\t" else ","
    read.delim(fs[1], sep = sep, check.names = FALSE,
               stringsAsFactors = FALSE)
  } else NULL
  if (!is.null(local)) names(local) <- make.unique(names(local))

  if (is.null(pr) && is.null(local)) return(NULL)
  if (is.null(local)) return(pr)
  if (is.null(pr))    return(local)

  local_fcol <- intersect(MSV_FILENAME_COLS, names(local))[1]
  if (is.na(local_fcol)) {
    message("Local TSV for ", deposit_id,
            " has no filename column; using Pan-ReDU only")
    return(pr)
  }

  # Case-insensitive stem join: sources differ in capitalization on identical stems.
  pr$.stem    <- tolower(tools::file_path_sans_ext(basename(pr$filename)))
  local$.stem <- tolower(tools::file_path_sans_ext(
                          basename(local[[local_fcol]])))

  add_cols <- setdiff(names(local), c(names(pr), ".stem"))
  out <- if (length(add_cols)) {
    merge(pr, local[, c(".stem", add_cols), drop = FALSE],
          by = ".stem", all.x = TRUE)
  } else pr

  # Submitter-only files, kept only when the deposit inventory confirms they exist.
  extra <- local[!local$.stem %in% pr$.stem, , drop = FALSE]
  if (nrow(extra)) {
    real <- tryCatch(
      tolower(tools::file_path_sans_ext(basename(massive_real_ms_files(deposit_id)))),
      error = function(e) {
        message("  massive_list_files failed for ", deposit_id,
                "; keeping Pan-ReDU rows only")
        NULL
      })
    extra <- if (is.null(real)) extra[0, , drop = FALSE]
             else extra[extra$.stem %in% real, , drop = FALSE]
  }
  if (nrow(extra)) {
    # Keep a submitter-only row only if >= MIN_SUBMITTER_FIELDS canonical fields resolve
    # to a real value (is_missing_val treats ReDU/MIxS placeholders as absent).
    n_fields <- rowSums(vapply(DRAFT_CANDIDATES, function(pat) {
      cl <- match_candidate(names(extra), pat)
      if (is.na(cl) || !nzchar(cl) || !cl %in% names(extra))
        rep(FALSE, nrow(extra)) else !is_missing_val(extra[[cl]])
    }, logical(nrow(extra))))
    keep  <- n_fields >= MIN_SUBMITTER_FIELDS
    thin  <- sum(!keep)
    extra <- extra[keep, , drop = FALSE]
    if (thin)
      message("  ", deposit_id, ": dropped ", thin,
              " submitter-only files with < ", MIN_SUBMITTER_FIELDS,
              " usable metadata fields")
  }
  if (nrow(extra)) {
    extra$filename <- extra[[local_fcol]]
    message("  ", deposit_id, ": +", nrow(extra),
            " submitter-only files not in Pan-ReDU")
    out <- dplyr::bind_rows(out, extra)
  }
  out$.stem <- NULL
  out
}

harvest_metadata_coverage <- function(deposit_id) {
  s <- native_sample_table(deposit_id)
  if (is.null(s) || !nrow(s)) {
    return(data.frame(
      dataset_id   = deposit_id,
      n_samples    = 0L,
      n_canonical_mapped = 0L,
      stringsAsFactors = FALSE
    ))
  }
  cols <- names(s)
  mapped <- vapply(names(DRAFT_CANDIDATES),
                   function(f) match_candidate(cols, DRAFT_CANDIDATES[[f]]),
                   character(1))

  # Demote a mapping to NA when its column holds only REDU/MIxS placeholders (see REDU_MISSING).
  has_value <- function(col_name) {
    if (is.na(col_name) || !col_name %in% names(s)) return(FALSE)
    any(!is_missing_val(s[[col_name]]))
  }
  mapped <- vapply(mapped,
                   function(c) if (has_value(c)) c else NA_character_,
                   character(1))

  row <- data.frame(
    dataset_id         = deposit_id,
    n_samples          = nrow(s),
    n_canonical_mapped = sum(!is.na(mapped)),
    stringsAsFactors   = FALSE
  )
  for (f in names(mapped)) row[[paste0("col_", f)]] <- unname(mapped[f])
  row
}

# ---- load -----------------------------------------------------------------

read_metadata_table <- function(
    path = "artifacts/dataset_metadata_table.csv") {
  if (!file.exists(path))
    stop("Metadata table not found: ", path,
         "\n  Render 01-masst-curation.qmd first; it builds this artifact.")
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Apply col_*/const_* spec -> one row per native sample, canonical cols (file joined per-assay later).
apply_spec_to_native <- function(native, spec) {
  out <- data.frame(row.names = seq_len(nrow(native)))
  for (canon in setdiff(CANONICAL_COLS, "file")) {
    col_key   <- paste0("col_",   canon)
    const_key <- paste0("const_", canon)
    const_val <- if (const_key %in% names(spec)) spec[[const_key]] else NA
    native_col <- if (col_key %in% names(spec)) spec[[col_key]] else NA

    if (!is.na(const_val) && nzchar(as.character(const_val))) {
      out[[canon]] <- const_val
    } else if (!is.na(native_col) && nzchar(native_col) &&
               native_col %in% names(native)) {
      v <- as.character(native[[native_col]])
      v[is_missing_val(v)] <- NA_character_
      out[[canon]] <- v
    } else {
      out[[canon]] <- NA
    }
  }
  out
}

# load_metadata(dataset_id): combined key; deposit_id and assay looked up from the table.
load_metadata <- function(dataset_id,
                          table_path =
                            "artifacts/dataset_metadata_table.csv") {
  tbl <- read_metadata_table(table_path)
  if (!dataset_id %in% tbl$dataset_id)
    stop("No row for ", dataset_id, " in ", table_path)
  spec       <- as.list(tbl[tbl$dataset_id == dataset_id, ][1, ])
  deposit_id <- spec$deposit_id
  assay      <- spec$assay

  native <- native_sample_table(deposit_id)
  if (is.null(native))
    stop("Could not fetch native sample metadata for ", deposit_id)
  s_canon <- apply_spec_to_native(native, spec)

  out <- if (startsWith(deposit_id, "MTBLS"))
           expand_files_mtbls(deposit_id, assay, native, s_canon, spec)
         else if (startsWith(deposit_id, "ST"))
           expand_files_mwb(deposit_id, assay, native, s_canon, spec)
         else if (startsWith(deposit_id, "MSV"))
           expand_files_massive(deposit_id, native, s_canon, spec)
         else
           stop("Unknown deposit_id prefix: ", deposit_id)

  for (k in CANONICAL_COLS) if (!k %in% names(out)) out[[k]] <- NA
  # `age` and `timepoint` are numeric by intent but come from free-text cells, and
  # deposited metadata puts non-numbers in them (MSV000082493 carries a literal
  # "FALSE" in ATTRIBUTE_Time_Point_Mins). Normalise here so every consumer sees a
  # number or NA, rather than each one re-deciding and a stray string reaching a CSV.
  for (k in c("age", "timepoint")) out[[k]] <- as_num(out[[k]])
  out[, CANONICAL_COLS]
}

expand_files_mtbls <- function(deposit_id, assay, native, s_canon, spec) {
  if (is.null(assay) || !nzchar(assay))
    stop("MTBLS load_metadata() requires `assay` (the a_*.txt file name).")
  # Prefer the copy `01` already cached in dataset_metadata/. Re-fetching here
  # made every render depend on EBI being reachable, and when its FTP timed out
  # the caller's tryCatch silently swapped a 35-column table for a 25-column
  # stub -- which surfaced two chunks later as an rbind column mismatch, not as
  # a network error. Downloads still go through the MsBackend package.
  cache <- file.path(META_DIR, assay)
  a_tbl <- if (file.exists(cache)) read.delim(cache, check.names = FALSE) else
    MsBackendMetaboLights::mtbls_assay_data(deposit_id, assayName = assay)
  fcol <- if ("Derived Spectral Data File" %in% names(a_tbl) &&
              any(nzchar(trimws(as.character(
                a_tbl[["Derived Spectral Data File"]])))))
            "Derived Spectral Data File"
          else
            "Raw Spectral Data File"
  files <- basename(trimws(as.character(a_tbl[[fcol]])))
  sample_name_col <- spec$col_sample_name
  if (is.na(sample_name_col) || !nzchar(sample_name_col))
    sample_name_col <- "Sample Name"
  per_file <- data.frame(
    file        = files,
    sample_name = a_tbl[[sample_name_col]],
    stringsAsFactors = FALSE
  )
  s_canon$sample_name <- native[[sample_name_col]]
  merge(per_file, s_canon[, c("sample_name", setdiff(names(s_canon),
                                                       "sample_name"))],
        by = "sample_name", all.x = TRUE, sort = FALSE)
}

expand_files_mwb <- function(deposit_id, assay, native, s_canon, spec) {
  if (is.null(assay) || !nzchar(assay))
    stop("ST load_metadata() requires `assay` (the AN file name).")
  md <- MsBackendMetabolomicsWorkbench::mwb_metadata(deposit_id)
  an_ids    <- md$MS_run$ANALYSIS_ID
  an_labels <- paste0(deposit_id, "_", an_ids, ".txt")
  an_idx    <- match(assay, an_labels)
  if (is.na(an_idx))
    stop("Assay ", assay, " not found in MS_run for ", deposit_id)
  rfn_col <- grep("RAW_FILE_NAME", names(native), value = TRUE)[1]
  if (is.na(rfn_col))
    stop("No RAW_FILE_NAME column in sample_annotation for ", deposit_id)
  per_sample <- strsplit(trimws(native[[rfn_col]]), "\\s+")
  files <- vapply(per_sample,
                  function(x) if (length(x) >= an_idx) x[an_idx] else NA,
                  character(1))
  out <- s_canon
  out$file <- files
  out$sample_name <- if (!is.na(spec$col_sample_name) &&
                         spec$col_sample_name %in% names(native))
                       native[[spec$col_sample_name]] else NA
  out[!is.na(out$file) & nzchar(out$file), ]
}

expand_files_massive <- function(deposit_id, native, s_canon, spec) {
  file_col <- match_candidate(names(native), MSV_FILENAME_COLS)
  if (is.na(file_col))
    stop("No filename column in MassIVE metadata for ", deposit_id)
  out <- s_canon
  out$file <- basename(as.character(native[[file_col]]))
  out$sample_name <- if (!is.na(spec$col_sample_name) &&
                         spec$col_sample_name %in% names(native))
                       native[[spec$col_sample_name]] else NA
  out
}

# ===== Pan-ReDU (was panredu.R) =====

# GNPS2's harmonised metadata layer over MassIVE/MetaboLights/Workbench; often the only
# standardised source for MassIVE deposits. Used as fallback/gap-fill alongside local TSVs.

PANREDU_URL_TMPL <- paste0(
  "https://redu.gnps2.org/attribute/ATTRIBUTE_DatasetAccession/",
  "attributeterm/%s/files?filters=%%5B%%5D"
)

# Fetch Pan-ReDU per-file metadata for one accession (MSV/MTBLS/ST). Returns a data.frame
# (one row per file) or NULL on failure/empty; caches to <cache_dir>/<ds>_panredu.tsv.
fetch_panredu_metadata <- function(ds,
                                   cache_dir = "dataset_metadata",
                                   force = FALSE,
                                   timeout = 180L) {
  cache_path <- file.path(cache_dir, paste0(ds, "_panredu.tsv"))
  if (!force && file.exists(cache_path)) {
    return(read.delim(cache_path, check.names = FALSE,
                      stringsAsFactors = FALSE, na.strings = ""))
  }

  url <- sprintf(PANREDU_URL_TMPL, ds)
  tmp <- tempfile(fileext = ".json")
  ok <- tryCatch({
    h <- curl::new_handle(ssl_verifypeer = 0, ssl_verifyhost = 0,
                          followlocation = 1, timeout = timeout)
    curl::curl_download(url, tmp, handle = h, quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Pan-ReDU fetch failed for ", ds, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    if (file.exists(tmp)) unlink(tmp)
    return(NULL)
  }

  parsed <- jsonlite::fromJSON(tmp, simplifyVector = TRUE)
  unlink(tmp)
  if ((is.list(parsed) && !length(parsed)) ||
      (is.data.frame(parsed) && !nrow(parsed))) {
    message("Pan-ReDU returned no rows for ", ds)
    return(NULL)
  }

  df <- as.data.frame(parsed, check.names = FALSE, stringsAsFactors = FALSE)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(df, cache_path, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "")
  message("Pan-ReDU cached -> ", cache_path, " (", nrow(df), " files, ",
          ncol(df), " columns)")
  df
}

# ---- extraction-set selection from the deposit inventory -----------------

# Map prefixed cache filenames back to the bare names the pipeline keys on. MsBackend prefixes
# cached files to dodge cross-dataset collisions: MetaboLights "<id>_<assayIndex>_<name>"
# (e.g. MTBLS1866_1_1_p.mzML), MassIVE "<id>_<name>", Workbench "<zipbase>_<name>". Detect the
# shared prefix from a file whose suffix is in `bare_ref`, then strip it; bare names pass through.
# Lives here (not ms1.R) because norm_file_key() calls it and scripts source metadata.R alone.
to_bare_name <- function(x, bare_ref) {
  if (!length(x)) return(x)
  bare_ref <- unique(bare_ref)
  # Resolve each name individually. Deriving one prefix from the first example assumes a
  # deposit has only one form, and MTBLS1866 carries both MTBLS1866_129_p.mzML and the
  # doubled MTBLS1866_1_MTBLS1866_DA44_p.mzML. Longest match wins, and the leading "_"
  # keeps 21_p from resolving against 1_p.
  vapply(x, function(f) {
    if (f %in% bare_ref) return(f)
    cand <- bare_ref[endsWith(f, paste0("_", bare_ref))]
    if (!length(cand)) return(f)
    cand[which.max(nchar(cand))]
  }, character(1), USE.NAMES = FALSE)
}

# Normalise a file name for cross-source matching: drop accession prefix + extension, lower-case,
# so MSV000082433_957.mzML, 957.mzXML and 957.mzML all reduce to "957".
# Pass `bare_ref` (the dataset's real bio_files) whenever `x` comes from the cache: MsBackend
# prefixes with accession and assay, and no regex can tell the assay token from the name, so
# to_bare_name() resolves it against the real file list. Omit it for repo-derived names.
#   cache-derived (dataOrigin, metrics_*$file) -> norm_file_key(x, bio_files)
#   repo-derived  (hits_curated$file_path, submitter TSVs) -> norm_file_key(x)
norm_file_key <- function(x, bare_ref = NULL) {
  b <- basename(as.character(x))
  if (!is.null(bare_ref) && length(bare_ref))
    b <- to_bare_name(b, basename(as.character(bare_ref)))
  b <- sub("^(MSV[0-9]+|MTBLS[0-9]+|ST[0-9]+)[_-]", "", b)
  tolower(sub("[.][^.]*$", "", b))
}

# Build the MS1 extraction set for a MassIVE deposit from the deposit's OWN inventory, not a
# metadata catalog: Pan-ReDU lists only harmonised-metadata files, a strict subset of what the
# deposit holds and what MASST searched (082433: 223 catalog vs 457 spectra).
# Returns data.frame(file, has_metadata); no-metadata files are kept only when they carry a
# confirmed hit and must be excluded from stratified analysis.
#   meta_keys = norm_file_key() of files with usable metadata (QC/blanks excluded by SampleType/EXCL_ST, not by name).
#   hit_files = basenames of confirmed-hit files; always kept, whatever tree they live in.
select_bio_files_massive <- function(deposit_id, meta_keys = character(),
                                     hit_files = character()) {
  inv <- tryCatch(MsBackendMassIVE::massive_list_files(deposit_id),
                  error = function(e) NULL)
  if (is.null(inv)) stop("massive_list_files failed for ", deposit_id)
  paths <- if (is.data.frame(inv)) {
    fc <- intersect(c("file_name", "fileName", "name", "path", "file",
                      "spectra_file"), names(inv))[1]
    as.character(inv[[fc]])
  } else as.character(inv)

  # One copy per sample, deduplicated by normalised name not by picking a single tree:
  # deposits hold the same sample under ccms_peak/, peak/, ...; restricting to one tree drops
  # collections added by later updates (082493's 341 blood_plasma files live only under
  # updates/.../peak/blood/, the CYP2C19 matrix). Prefer the ccms_peak copy on duplicates.
  ms  <- grep("[.](mzML|mzXML)$", paths, value = TRUE, ignore.case = TRUE)
  ms  <- ms[order(!grepl("^ccms_peak/", ms, ignore.case = TRUE))]
  bio <- ms[!duplicated(norm_file_key(basename(ms)))]
  with_meta <- basename(bio)[norm_file_key(basename(bio)) %in% meta_keys]
  # Plus every confirmed-hit file, whatever tree it sits in.
  all_files <- unique(c(with_meta, basename(as.character(hit_files))))

  # Pooled QC dropped outright (EXCL_POOL): not an independent sample, so out of num and denom.
  n_pool <- sum(grepl(EXCL_POOL, all_files, perl = TRUE))
  if (n_pool) {
    message("  ", deposit_id, ": dropped ", n_pool, " pooled-QC files")
    all_files <- all_files[!grepl(EXCL_POOL, all_files, perl = TRUE)]
  }

  data.frame(file         = all_files,
             has_metadata = norm_file_key(all_files) %in% meta_keys,
             stringsAsFactors = FALSE)
}

# ===== assay/technical metadata (was assay_metadata.R) =====

# fetch_assay_metadata(dataset_ids): one row per (dataset_id, assay). Sources: MTBLS from
# mtbls_assay_data, ST from mwb_metadata()$MS_run, MSV from local TSV then Pan-ReDU fallback.

# ReDU/MIxS placeholders treated as missing.
REDU_MISSING <- c("missing value", "not applicable", "not collected",
                  "not provided", "restricted access", "na", "n/a", "",
                  "-", "unknown", "no data")

# Numeric where the text is a number, NA where it is not.
#
# as.numeric() answers "is this numeric?" by raising a warning, which then has to
# be hidden -- and hiding it is how a corrupt cell stays invisible. Deposited
# metadata is full of non-numeric entries in numeric columns (MSV000082493 carries
# a literal "FALSE" in ATTRIBUTE_Time_Point_Mins), so decide with a predicate and
# convert only what passes.
as_num <- function(x) {
  x <- as.character(x)
  ok <- !is.na(x) & grepl("^ *[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)? *$", x)
  out <- rep(NA_real_, length(x))
  out[ok] <- as.numeric(x[ok])
  out
}

is_missing_val <- function(x) {
  is.na(x) | tolower(trimws(as.character(x))) %in% REDU_MISSING
}

chrom_type_of <- function(...) {
  x <- paste(unlist(list(...)), collapse = " | ")
  if (grepl("HILIC|ZIC|pHILIC|hydrophilic", x, ignore.case = TRUE)) "HILIC"
  else if (grepl("\\bGC\\b|gas chrom",      x, ignore.case = TRUE)) "GC"
  else if (grepl("RP|C18|reverse",          x, ignore.case = TRUE)) "RP"
  else                                                              "Unknown"
}

first_nonempty <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nchar(x) > 0]
  if (length(x)) x[1] else NA_character_
}

fetch_assay_metadata_mtbls <- function(ds) {
  files <- MsBackendMetaboLights::mtbls_list_files(ds)
  assay_files <- grep("^a_.*\\.txt$", basename(as.character(files)),
                      value = TRUE)
  if (!length(assay_files)) return(NULL)

  do.call(rbind, lapply(assay_files, function(an) {
    ad <- MsBackendMetaboLights::mtbls_assay_data(ds, assayName = an)
    if (is.null(ad)) return(NULL)
    col_type  <- first_nonempty(ad[["Parameter Value[Column type]"]])
    col_model <- first_nonempty(ad[["Parameter Value[Column model]"]])
    chr_inst  <- first_nonempty(
      ad[["Parameter Value[Chromatography Instrument]"]])
    data.frame(
      dataset_id    = ds,
      assay         = an,
      repository    = "MetaboLights",
      chrom_type    = chrom_type_of(col_type, col_model, chr_inst, an),
      column_model  = col_model,
      polarity      = tolower(first_nonempty(
        ad[["Parameter Value[Scan polarity]"]])),
      mz_range      = first_nonempty(ad[["Parameter Value[Scan m/z range]"]]),
      instrument    = first_nonempty(ad[["Parameter Value[Instrument]"]]),
      mass_analyzer = first_nonempty(ad[["Parameter Value[Mass analyzer]"]]),
      ion_source    = first_nonempty(ad[["Parameter Value[Ion source]"]]),
      stringsAsFactors = FALSE
    )
  }))
}

fetch_assay_metadata_mwb <- function(ds) {
  md <- MsBackendMetabolomicsWorkbench::mwb_metadata(ds)
  if (is.null(md) || is.null(md$MS_run) || !nrow(md$MS_run)) return(NULL)
  m <- md$MS_run
  data.frame(
    dataset_id    = ds,
    assay         = paste0(ds, "_", m$ANALYSIS_ID, ".txt"),
    repository    = "Metabolomics Workbench",
    chrom_type    = vapply(m$CHROMATOGRAPHY_TYPE, chrom_type_of, character(1)),
    column_model  = m$COLUMN_NAME,
    polarity      = tolower(m$ION_MODE),
    mz_range      = NA_character_,
    instrument    = m$INSTRUMENT_NAME,
    mass_analyzer = m$INSTRUMENT_TYPE,
    ion_source    = m$MS_TYPE,
    stringsAsFactors = FALSE
  )
}

fetch_assay_metadata_massive <- function(ds) {
  clean <- function(x) {
    if (is.null(x)) return(NULL)
    ifelse(is_missing_val(x), NA_character_, as.character(x))
  }
  get_col <- function(df, col) {
    if (is.null(df) || !col %in% names(df)) NULL else clean(df[[col]])
  }
  all_missing <- function(x) is.null(x) || all(is.na(x))

  fs <- list.files("dataset_metadata",
                   pattern = paste0("^", ds, "_metadata\\.(csv|tsv)$"),
                   full.names = TRUE)
  local <- if (length(fs)) {
    sep <- if (endsWith(fs[1], ".tsv")) "\t" else ","
    read.delim(fs[1], sep = sep, check.names = FALSE,
               stringsAsFactors = FALSE)
  } else NULL

  chr <- get_col(local, "ChromatographyAndPhase")
  ion <- get_col(local, "IonizationSourceAndPolarity")
  ms  <- get_col(local, "MassSpectrometer")

  if (all_missing(chr) || all_missing(ion) || all_missing(ms)) {
    pr <- fetch_panredu_metadata(ds)
    if (!is.null(pr) && nrow(pr)) {
      if (all_missing(chr)) chr <- get_col(pr, "ChromatographyAndPhase")
      if (all_missing(ion)) ion <- get_col(pr, "IonizationSourceAndPolarity")
      if (all_missing(ms))  ms  <- get_col(pr, "MassSpectrometer")
    }
  }

  if (all_missing(chr) && all_missing(ion) && all_missing(ms)) {
    return(data.frame(
      dataset_id    = ds, assay = ds, repository = "MassIVE",
      chrom_type    = "Unknown",
      column_model  = NA_character_, polarity = NA_character_,
      mz_range      = NA_character_, instrument = NA_character_,
      mass_analyzer = NA_character_, ion_source = NA_character_,
      stringsAsFactors = FALSE))
  }

  # one row per unique (chrom, ion, instrument) combo
  combos <- unique(data.frame(chr = chr, ion = ion, ms = ms,
                              stringsAsFactors = FALSE))
  do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
    data.frame(
      dataset_id    = ds,
      assay         = paste0(ds, "_combo", i),
      repository    = "MassIVE",
      chrom_type    = chrom_type_of(combos$chr[i]),
      column_model  = combos$chr[i],
      polarity      = if (!is.na(combos$ion[i]))
                        tolower(sub(".*\\(([^)]+)\\).*", "\\1", combos$ion[i]))
                      else NA_character_,
      mz_range      = NA_character_,
      instrument    = if (!is.na(combos$ms[i]))
                        sub("\\|.*", "", combos$ms[i]) else NA_character_,
      mass_analyzer = NA_character_,
      ion_source    = if (!is.na(combos$ion[i]))
                        sub(" \\(.*", "", combos$ion[i]) else NA_character_,
      stringsAsFactors = FALSE)
  }))
}

fetch_assay_metadata <- function(dataset_ids) {
  rows <- lapply(unique(dataset_ids), function(ds) {
    if      (startsWith(ds, "MSV"))   fetch_assay_metadata_massive(ds)
    else if (startsWith(ds, "MTBLS")) fetch_assay_metadata_mtbls(ds)
    else if (startsWith(ds, "ST"))    fetch_assay_metadata_mwb(ds)
    else NULL
  })
  do.call(rbind, rows)
}

# ===== biological file selection (was bio_files.R) =====

META_DIR <- "dataset_metadata"
EXCL_ST  <- "(?i)^(blank|qc|pool|reference|standard|control|empty)"

# Pooled QC excluded by file name: 094097's "Pool_*" files are SampleType "animal" so EXCL_ST
# never fires, yet a re-injected pool inflates counts and denominator. The one name-based
# override; deliberately narrow (leading "pool" only).
EXCL_POOL <- "(?i)^pool"

# Assay file names for a MetaboLights study from its ISA investigation file (cached locally).
mtbls_inv_assay_names <- function(ds) {
  inv_path <- file.path(META_DIR, paste0(ds, "_i_Investigation.txt"))
  if (!file.exists(inv_path)) {
    inv_url   <- paste0("https://ftp.ebi.ac.uk/pub/databases/metabolights/",
                        "studies/public/", ds, "/i_Investigation.txt")
    inv_lines <- readLines(url(inv_url), warn = FALSE)
    writeLines(inv_lines, inv_path)
  } else {
    inv_lines <- readLines(inv_path, warn = FALSE)
  }
  assay_line <- inv_lines[startsWith(inv_lines, "Study Assay File Name")]
  if (!length(assay_line)) return(NULL)
  an <- trimws(unlist(strsplit(sub("^[^\t]+\t", "", assay_line[1]), "\t")))
  an[nchar(an) > 0]
}

# Fetch (or load from cache) one MetaboLights assay table; tries mtbls_assay_data() then FTP.
fetch_assay_table <- function(ds, assay_name) {
  cache_path <- file.path(META_DIR, assay_name)
  if (file.exists(cache_path))
    return(read.delim(cache_path, check.names = FALSE))
  adf <- tryCatch(mtbls_assay_data(ds, assay_name),
                  error = function(e) {
                    message("  mtbls_assay_data failed for ", ds, "/",
                            assay_name, " (", e$message,
                            "); falling back to FTP")
                    NULL
                  })
  if (!is.null(adf)) {
    write.table(adf, cache_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  Assay cached (mtbls_assay_data) → ", cache_path)
    return(adf)
  }
  assay_url <- paste0("https://ftp.ebi.ac.uk/pub/databases/metabolights/",
                      "studies/public/", ds, "/", assay_name)
  adf <- read.delim(url(assay_url), check.names = FALSE)
  write.table(adf, cache_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  Assay cached (FTP) → ", cache_path)
  adf
}

# Pick the spectral file column from an ISA assay table: prefer Derived over Raw (Raw is often
# all-NA while Derived has paths). NA if neither exists (e.g. NMR-only assay).
fcol_pick <- function(adf) {
  derived <- grep("Derived Spectral Data File", names(adf),
                  value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(derived) && any(!is.na(adf[[derived]]))) return(derived)
  grep("Raw Spectral Data File", names(adf),
       value = TRUE, ignore.case = TRUE)[1]
}

# Named list assay_name -> file basenames; NULL for MassIVE (no assay concept), NMR assays skipped.
assay_lookup_for <- function(ds) {
  if (startsWith(ds, "MTBLS")) {
    an_names <- mtbls_inv_assay_names(ds)
    if (is.null(an_names) || !length(an_names)) return(NULL)
    out <- list()
    for (an in an_names) {
      adf <- fetch_assay_table(ds, an)
      if (is.null(adf) || !nrow(adf)) next
      fc <- fcol_pick(adf)
      if (is.na(fc) || !fc %in% names(adf)) next     # NMR or missing column
      files <- basename(trimws(as.character(adf[[fc]])))
      files <- files[!is.na(files) & nzchar(files)]
      if (length(files)) out[[an]] <- unique(files)
    }
    if (length(out)) out else NULL
  } else if (startsWith(ds, "ST")) {
    md <- MsBackendMetabolomicsWorkbench::mwb_metadata(ds)
    if (is.null(md) || is.null(md$MS_run) || !nrow(md$MS_run)) return(NULL)
    an_ids    <- md$MS_run$ANALYSIS_ID
    an_labels <- paste0(ds, "_", an_ids, ".txt")
    rfn_col   <- grep("RAW_FILE_NAME", names(md$sample_annotation),
                      value = TRUE)[1]
    if (is.na(rfn_col)) return(NULL)
    per_sample <- strsplit(trimws(md$sample_annotation[[rfn_col]]), "\\s+")
    out <- list()
    for (i in seq_along(an_labels)) {
      files <- vapply(per_sample,
                      function(x) if (length(x) >= i) x[i] else NA_character_,
                      character(1))
      files <- basename(files[!is.na(files) & nzchar(files)])
      if (length(files)) out[[an_labels[i]]] <- unique(files)
    }
    if (length(out)) out else NULL
  } else NULL
}

# Which assay a file belongs to in a lookup; NA if no assay's list has this basename.
file_to_assay <- function(file_basename, lookup) {
  if (is.null(lookup) || is.na(file_basename) || !nzchar(file_basename))
    return(NA_character_)
  for (an in names(lookup)) {
    if (file_basename %in% lookup[[an]]) return(an)
  }
  NA_character_
}

# MassIVE phantom-file guard: Pan-ReDU can advertise .mzXML names never deposited as spectra;
# such phantoms made massive_sync_data_files() abort the dataset (084008: 21 phantoms). Keep only
# files truly present as .mzML/.mzXML in the deposit listing. MsBackendMassIVE only (never curl/FTP).
massive_real_ms_files <- function(deposit_id) {
  fl <- MsBackendMassIVE::massive_list_files(deposit_id)
  nm <- if (is.data.frame(fl)) {
    unlist(fl[[grep("file|name|path", names(fl), ignore.case = TRUE)[1]]])
  } else fl
  bn <- basename(as.character(nm))
  bn[grepl("[.]mz(X)?ML$", bn, ignore.case = TRUE)]
}
