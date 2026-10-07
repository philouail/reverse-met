# Per-dataset biological metadata: harvest_metadata_coverage() matches canonical to native
# columns; load_metadata() applies the table's col_*/const_* spec to the native table.

CANONICAL_COLS <- c("file", "sample_name", "subject_id", "sex", "age",
                    "group", "health_status", "body_site", "country",
                    "treatment", "timepoint")

# MIN_SUBMITTER_FIELDS: canonical fields a submitter-only file must fill to be kept.
# MSV_FILENAME_COLS: native raw-file-name columns (naming varies across submitter TSVs).
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

# Display label for a dataset_id: the accession without the assay boilerplate;
# ids that would share a label keep their full id.
short_ds <- function(x, ids = unique(x)) {
  ids <- unique(as.character(ids))
  s   <- sub("_AN[0-9]+$", "", sub("_LC-MS_.*$", "", ids))
  lab <- stats::setNames(ifelse(s %in% s[duplicated(s)], ids, s), ids)
  out <- unname(lab[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

# Candidate native columns per canonical field, matched case-insensitively on the full
# name, first hit wins (ISA-Tab, mwTab, MassIVE ATTRIBUTE_*, Pan-ReDU).
DRAFT_CANDIDATES <- list(
  subject_id    = c("Source Name", "Subject", "Subject ID", "SubjectID",
                    "ATTRIBUTE_Subject",
                    "ATTRIBUTE_SubjectIdentifierAsRecorded",
                    "SubjectIdentifierAsRecorded", "UniqueSubjectID",
                    ".*_SubjectIdentifierAsRecorded",
                    "ATTRIBUTE_Subject.*", "subject.?id", "subject",
                    "patient.?id", "participant.?id", "donor.?id"),
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
                    "UBERONBodyPartName", ".*BodyPartOntologyName",
                    "ATTRIBUTE_Tissue", "ATTRIBUTE_Matrix", "ATTRIBUTE_Biofluid",
                    "tissue", "matrix", "biofluid", "body.?site", "sample_type",
                    "specimen.?type", "specimen"),
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
                    "ATTRIBUTE_Sex", "ATTRIBUTE_Subject_Sex",
                    "BiologicalSex",
                    "ATTRIBUTE_Gender", "sex", "gender"),
  age           = c("Characteristics\\[Age\\]",
                    "Factor Value\\[Age\\]",
                    "Comment\\[Age\\]",
                    "Additional sample data: Age",
                    "ATTRIBUTE_Age",
                    "AgeInYears", ".*_AgeInYears",
                    "ATTRIBUTE_Age.*", "age", "age.?\\(.*\\)", "age.?years?",
                    "age_.*", ".*_age"),
  country       = c("Comment\\[Country\\]",
                    "Characteristics\\[Geographic location\\]",
                    "Characteristics\\[Country\\]",
                    "ATTRIBUTE_Country",
                    "Country"),
  # Species: Pan-ReDU, ISA-Tab, the Workbench study (SUBJECT_SPECIES) and
  # submitter tables. Not a canonical column: step 2 reads it for its inventory.
  species       = c("NCBITaxonomy", "Characteristics\\[Organism\\]",
                    "SUBJECT_SPECIES", "ATTRIBUTE_Organism", "ATTRIBUTE_Species",
                    "organism", "species"),
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

# ---- repositories that do not answer --------------------------------------

# Error-message words of a repository that did not answer (connection failure, timeout,
# HTTP 429/502/503/504); a 404 is not one. Keep in sync with NETWORK in bench/rerun.sh.
NETWORK_ERROR <- paste(
  "Failed to perform HTTP request", "resolve host", "connect to server",
  "Timeout (of [0-9]+ seconds )?was reached", "timed out", "Connection (was )?reset",
  "Connection refused", "Recv failure", "Send failure", "receiving data",
  "sending data", "Server returned nothing", "Transferred a partial file",
  "SSL connect error", "(HTTP|error:?) ?(429|50[234])", "Too Many Requests",
  "Bad Gateway", "Service Unavailable", "Gateway Time-?out", "not all .rnames. found",
  sep = "|")

# TRUE when `e` (error or warning) is a repository that did not answer; handlers re-raise
# it so the render stops (and bench/rerun.sh retries) instead of writing a smaller result.
is_network_error <- function(e)
  grepl(NETWORK_ERROR, conditionMessage(e), ignore.case = TRUE)

# url() (also behind read.table() on a URL) reports a network failure only as a warning;
# raise that warning as the error.
offline_as_error <- function(expr)
  withCallingHandlers(expr, warning = function(w)
    if (is_network_error(w)) stop(conditionMessage(w), call. = FALSE))

# A study deposited twice (same_study_different_arm) counts each file once: TRUE for files
# of `deposits` also in the relation's first-column deposit; files_of(d) lists d's files.
is_copy <- function(files, deposits, files_of,
                    related_csv = "artifacts/related_deposits.csv") {
  out <- logical(length(files))
  if (!file.exists(related_csv)) return(out)
  r <- read.csv(related_csv, stringsAsFactors = FALSE)
  r <- r[r$relation == "same_study_different_arm", , drop = FALSE]
  for (i in seq_len(nrow(r))) {
    in_copy <- deposits == r$related_to[i]
    out[in_copy] <- norm_file_key(files[in_copy]) %in% norm_file_key(files_of(r$dataset_id[i]))
  }
  out
}

# ---- harvest --------------------------------------------------------------

native_sample_table <- function(deposit_id) {
  if (startsWith(deposit_id, "MTBLS")) {
    # Cached ISA-Tab study table first, so a render does not depend on EBI being reachable;
    # step 5 reads the s_<id>.txt copy too.
    cache <- file.path(META_DIR, paste0("s_", deposit_id, ".txt"))
    if (!file.exists(cache)) {
      s <- offline_as_error(MsBackendMetaboLights::mtbls_sample_data(deposit_id))
      dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)
      write.table(s, cache, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    read.delim(cache, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (startsWith(deposit_id, "ST")) {
    # The cached copy first, as for MetaboLights (mwb_metadata_cached()).
    md <- mwb_metadata_cached(deposit_id)
    s  <- md$sample_annotation
    # Species and sample type are stated once per study, not per sample.
    s$SUBJECT_SPECIES <- md$MS_run$SUBJECT_SPECIES[1]
    s$SAMPLE_TYPE     <- md$MS_run$SAMPLE_TYPE[1]
    s
  } else if (startsWith(deposit_id, "MSV")) {
    msv_native_sample_table(deposit_id)
  } else NULL
}

# MassIVE: Pan-ReDU joined with the submitter TSV on file stem; submitter-only rows need the
# file in the deposit listing and min_fields filled (0 keeps all). NULL if neither exists.
msv_native_sample_table <- function(deposit_id, min_fields = MIN_SUBMITTER_FIELDS) {
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

  if (is.null(pr) && is.null(local)) {
    # Nothing found is not the same as Pan-ReDU not answering.
    if (panredu_failed(deposit_id)) stop("Pan-ReDU did not answer")
    return(NULL)
  }
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
        if (is_network_error(e)) stop(e)
        message("  massive_list_files failed for ", deposit_id,
                "; keeping Pan-ReDU rows only")
        NULL
      })
    extra <- if (is.null(real)) extra[0, , drop = FALSE]
             else extra[extra$.stem %in% real, , drop = FALSE]
  }
  if (nrow(extra)) {
    # Keep a submitter-only row only if >= min_fields canonical fields resolve
    # to a real value (is_missing_val treats ReDU/MIxS placeholders as absent).
    n_fields <- rowSums(vapply(DRAFT_CANDIDATES, function(pat) {
      cl <- match_candidate(names(extra), pat)
      if (is.na(cl) || !nzchar(cl) || !cl %in% names(extra))
        rep(FALSE, nrow(extra)) else !is_missing_val(extra[[cl]])
    }, logical(nrow(extra))))
    keep  <- n_fields >= min_fields
    thin  <- sum(!keep)
    extra <- extra[keep, , drop = FALSE]
    if (thin)
      message("  ", deposit_id, ": dropped ", thin,
              " submitter-only files with < ", min_fields,
              " usable metadata fields")
  }
  if (nrow(extra)) {
    extra$filename <- extra[[local_fcol]]
    message("  ", deposit_id, ": +", nrow(extra),
            " submitter-only files not in Pan-ReDU")
    # The two sources can type the same column differently (a year as text in one,
    # a number in the other), which bind_rows refuses; stack them as text.
    as_text <- function(d) { d[] <- lapply(d, as.character); d }
    out <- dplyr::bind_rows(as_text(out), as_text(extra))
  }
  out$.stem <- NULL
  out
}

harvest_metadata_coverage <- function(deposit_id) {
  # Unreadable metadata gives n_samples NA (not 0) so step 3 can tell it from "nothing
  # found"; a repository that does not answer stops the render.
  s <- tryCatch(native_sample_table(deposit_id), error = function(e) {
    if (is_network_error(e)) stop(e)
    message("  metadata not reachable for ", deposit_id, ": ", conditionMessage(e))
    FALSE
  })
  if (isFALSE(s))
    return(data.frame(dataset_id = deposit_id, n_samples = NA_integer_,
                      n_canonical_mapped = 0L, stringsAsFactors = FALSE))
  if (is.null(s) || !nrow(s)) {
    return(data.frame(
      dataset_id   = deposit_id,
      n_samples    = 0L,
      n_canonical_mapped = 0L,
      stringsAsFactors = FALSE
    ))
  }
  # For each field, the first candidate column that holds a real value: a column of
  # REDU/MIxS placeholders only (see REDU_MISSING) passes to the next candidate.
  has_value <- function(col_name) any(!is_missing_val(s[[col_name]]))
  mapped <- vapply(names(DRAFT_CANDIDATES), function(f) {
    for (pat in DRAFT_CANDIDATES[[f]]) {
      hits <- grep(paste0("^", pat, "$"), names(s), value = TRUE, ignore.case = TRUE)
      hits <- hits[vapply(hits, has_value, logical(1))]
      if (length(hits)) return(hits[1])
    }
    NA_character_
  }, character(1))

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
         "\n  Render 03-curation.qmd first; it builds this artifact.")
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
  # A deposit that suffixes each subject id with the visit (MSV000084556: FD001_T1,
  # FD001_T2, ...) declares the suffix in step 2, so that one person is one subject.
  if (!is.null(spec$const_subject_suffix) && !is.na(spec$const_subject_suffix))
    out$subject_id <- sub(spec$const_subject_suffix, "", out$subject_id)
  # A deposit that names each person only in its file names declares in step 2 a pattern
  # whose first bracket is the person (MSV000094097: OAD_001_Basal_P).
  pat <- spec$const_subject_from_file
  if (!is.null(pat) && !is.na(pat)) {
    fcol <- match_candidate(names(native), MSV_FILENAME_COLS)
    if (!is.na(fcol)) {
      f   <- basename(as.character(native[[fcol]]))
      hit <- is.na(out$subject_id) & grepl(pat, f)
      out$subject_id[hit] <- sub(pat, "\\1", f[hit])
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
  # age and timepoint come from free-text cells that can hold non-numbers (a literal "FALSE"):
  # normalise here so every consumer sees a number or NA.
  for (k in c("age", "timepoint")) out[[k]] <- as_num(out[[k]])
  out[, CANONICAL_COLS]
}

expand_files_mtbls <- function(deposit_id, assay, native, s_canon, spec) {
  if (is.null(assay) || !nzchar(assay))
    stop("MTBLS load_metadata() requires `assay` (the a_*.txt file name).")
  # Prefer the copy `02` already cached in dataset_metadata/, so that a render does not
  # depend on EBI being reachable; otherwise the table is fetched with MsBackendMetaboLights.
  cache <- file.path(META_DIR, assay)
  a_tbl <- if (file.exists(cache)) read.delim(cache, check.names = FALSE) else
    offline_as_error(MsBackendMetaboLights::mtbls_assay_data(deposit_id, assayName = assay))
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
  an_labels <- paste0(deposit_id, "_", mwb_analysis_order(deposit_id), ".txt")
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

# ===== what ReDU already knows (step 2, first pass) =====

# Deposits removed on ReDU alone, judged only when ReDU describes every hit file (a missing species
# or body site never removes one): not_human, no_blood (with `body_site`), not_reachable, `by_hand`.
redu_first_pass <- function(all_hits, fasst_db, redu_out, body_site = NULL,
                            by_hand = NULL) {
  deposits <- unique(all_hits$ATTRIBUTE_DatasetAccession)
  if (file.exists(fasst_db)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), fasst_db, flags = RSQLite::SQLITE_RO)
    r <- DBI::dbGetQuery(con, "SELECT ATTRIBUTE_DatasetAccession AS deposit_id, filename,
                                      NCBITaxonomy AS species, UBERONBodyPartName AS body_site,
                                      SampleType AS sample_type FROM redu_table")
    DBI::dbDisconnect(con)
    r$deposit_id <- trimws(gsub('"', "", r$deposit_id))
    write.csv(r[r$deposit_id %in% deposits, ], redu_out, row.names = FALSE)
  }
  r <- read.csv(redu_out, stringsAsFactors = FALSE)
  r <- r[r$deposit_id %in% deposits, ]

  # Deposits with a hit file that ReDU does not describe: not judged here.
  file_key  <- function(f) tolower(sub("[.][^.]*$", "", basename(f)))
  hit_files <- unique(all_hits[, c("ATTRIBUTE_DatasetAccession", "file_path")])
  described <- paste(hit_files$ATTRIBUTE_DatasetAccession, file_key(hit_files$file_path)) %in%
               paste(r$deposit_id, file_key(r$filename))
  partly    <- unique(hit_files$ATTRIBUTE_DatasetAccession[!described])

  r <- r[!grepl("blank|qc|standard", r$sample_type, ignore.case = TRUE), ]
  n_described <- length(setdiff(unique(r$deposit_id), partly))
  known <- data.frame(deposit_id = character(), reason = character())
  if (nrow(r)) {
    r$described     <- 1
    r$human         <- grepl("9606|homo sapiens", r$species, ignore.case = TRUE)
    r$other_species <- !is_missing_val(r$species) & !r$human
    # Without `body_site` no file counts as another site, so no deposit is no_blood.
    if (is.null(body_site)) {
      r$blood <- FALSE
      r$other_site <- FALSE
    } else {
      r$blood      <- grepl(body_site, r$body_site, ignore.case = TRUE)
      r$other_site <- !is_missing_val(r$body_site) & !r$blood
    }
    known <- aggregate(cbind(described, human, other_species, blood, other_site) ~ deposit_id,
                       data = r, FUN = sum)
    known <- known[!known$deposit_id %in% partly, ]
    known$reason <- ifelse(known$other_species > 0 & known$human == 0, "not_human",
                    ifelse(known$blood == 0 & known$other_site == known$described,
                           "no_blood", NA_character_))
  }
  unreadable <- deposits[!grepl("^(MSV|MTBLS|ST)", deposits)]
  if (!is.null(by_hand)) by_hand <- by_hand[by_hand$dataset_id %in% deposits, ]
  removed <- dplyr::bind_rows(
    tibble::tibble(dataset_id = unreadable, reason = "not_reachable",
                   notes      = "a repository the MsBackend packages do not read"),
    tibble::tibble(dataset_id = known$deposit_id[!is.na(known$reason)],
                   reason     = known$reason[!is.na(known$reason)],
                   notes      = "ReDU describes every file with a hit"),
    by_hand)
  removed <- removed[!duplicated(removed$dataset_id), ]
  left    <- setdiff(deposits, removed$dataset_id)
  list(removed = removed, deposits = left, n_described = n_described,
       hits = all_hits[all_hits$ATTRIBUTE_DatasetAccession %in% left, ])
}

# ===== Pan-ReDU =====

# GNPS2's harmonised metadata layer over MassIVE/MetaboLights/Workbench; often the only
# standardised source for MassIVE deposits. Used as fallback/gap-fill alongside local TSVs.

PANREDU_URL_TMPL <- paste0(
  "https://redu.gnps2.org/attribute/ATTRIBUTE_DatasetAccession/",
  "attributeterm/%s/files?filters=%%5B%%5D"
)

# Accessions Pan-ReDU failed for in this R session (an HTTP error, not a network failure),
# so each is asked once per render.
PANREDU_FAILED <- new.env()
panredu_failed <- function(ds) exists(ds, envir = PANREDU_FAILED, inherits = FALSE)

# Pan-ReDU per-file metadata for one accession, or NULL on failure/no rows; cached as
# <cache_dir>/<ds>_panredu.tsv, or an empty <ds>_panredu.none when Pan-ReDU has no rows.
fetch_panredu_metadata <- function(ds,
                                   cache_dir = "dataset_metadata",
                                   force = FALSE,
                                   timeout = 180L) {
  cache_path <- file.path(cache_dir, paste0(ds, "_panredu.tsv"))
  none_path  <- file.path(cache_dir, paste0(ds, "_panredu.none"))
  if (!force && file.exists(cache_path)) {
    return(read.delim(cache_path, check.names = FALSE,
                      stringsAsFactors = FALSE, na.strings = ""))
  }
  if (!force && (file.exists(none_path) || panredu_failed(ds))) return(NULL)

  url <- sprintf(PANREDU_URL_TMPL, ds)
  tmp <- tempfile(fileext = ".json")
  ok <- tryCatch({
    h <- curl::new_handle(ssl_verifypeer = 0, ssl_verifyhost = 0,
                          followlocation = 1, timeout = timeout)
    curl::curl_download(url, tmp, handle = h, quiet = TRUE)
    TRUE
  }, error = function(e) {
    if (is_network_error(e)) stop(e)
    message("Pan-ReDU fetch failed for ", ds, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    if (file.exists(tmp)) unlink(tmp)
    assign(ds, TRUE, envir = PANREDU_FAILED)
    return(NULL)
  }

  parsed <- jsonlite::fromJSON(tmp, simplifyVector = TRUE)
  unlink(tmp)
  if ((is.list(parsed) && !length(parsed)) ||
      (is.data.frame(parsed) && !nrow(parsed))) {
    message("Pan-ReDU returned no rows for ", ds)
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    file.create(none_path)
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

# MS data files of a MassIVE deposit, by extension; vendor folders (x.d/, x.raw/)
# count as one file each.
MSV_MS_EXT <- "[.](mzml|mzxml|mgf|raw|wiff|cdf|mzdata)$"
MSV_VENDOR_DIR <- "^(.*/)?([^/]+)[.](d|raw)/.*$"

# Candidate tables downloaded per deposit, at most: a few deposits carry hundreds of
# tables, and the ranking puts a sample table near the top.
MSV_MAX_TABLES <- 20L

# Tables that are not sample tables, by folder and by file name. Generic words like "method"
# match the file name only, since a methods/ folder can hold the sample table.
MSV_NOT_SAMPLE_DIR <- paste(
  "^search", "sirius", "csi", "canopus", "zodiac", "trinity", "proteosafe", "ili_files",
  "pdata", "quant", "feature", "mzmine", "msdial", "xcms", "standard",
  sep = "|")
MSV_NOT_SAMPLE_NAME <- paste(
  "quant", "feature", "edges", "cluster", "network", "sirius", "csi", "canopus",
  "zodiac", "formula", "structure", "identification", "trinity", "mzmine", "msdial",
  "xcms", "peaklist", "peaktable", "_ili", "lcparms", "nugenesis", "acqmethod",
  "method", "standard", "std_?mix", "compound", "smiles", "mwtab", "_maf", "read ?me",
  "msgs", "result", "audit",
  sep = "|")

# A file's key: its base name without extension, in lower case ("C:\\run\\A1.mzML"
# and "a1" are the same file).
msv_file_key <- function(x) {
  tolower(tools::file_path_sans_ext(sub(".*[/\\\\]", "", trimws(x))))
}

# Keys of the deposit's MS files, from its file listing.
msv_ms_file_keys <- function(fl) {
  in_vendor <- grepl(MSV_VENDOR_DIR, fl, ignore.case = TRUE)
  unique(c(msv_file_key(fl[!in_vendor & grepl(MSV_MS_EXT, fl, ignore.case = TRUE)]),
           tolower(sub(MSV_VENDOR_DIR, "\\2", fl[in_vendor], ignore.case = TRUE))))
}

# Table files that may be the sample table, most likely first: "meta" names, then sample-like
# names, then the rest, later updates first; a name with "metadata" is never dropped.
msv_table_candidates <- function(fl) {
  cand <- grep("[.](tsv|txt|csv|xlsx|xls)$", fl, value = TRUE, ignore.case = TRUE)
  cand <- cand[!grepl("[.](d|raw)/|(^|/)ccms_", cand, ignore.case = TRUE)]
  rel  <- sub("^updates/[^/]+/", "", tolower(cand))
  nm   <- basename(rel)
  drop <- grepl(MSV_NOT_SAMPLE_DIR, dirname(rel)) |
          (grepl(MSV_NOT_SAMPLE_NAME, nm) & !grepl("metadata", nm))
  cand <- cand[!drop]
  nm   <- nm[!drop]
  tier <- ifelse(grepl("meta(?!bol)", nm, perl = TRUE), 1,
          ifelse(grepl("sample|design|mapping|attribute|annot|group|info|key|manifest",
                       nm), 2, 3))
  upd  <- ifelse(startsWith(cand, "updates/"), sub("^updates/([^/]+)/.*$", "\\1", cand), "")
  head(cand[order(tier, -rank(upd))], MSV_MAX_TABLES)
}

# A downloaded table (first Excel sheet, or tab/comma/semicolon text) as text columns, or NULL.
# The header is the first of the top ten rows with >= 2 cells and > half the median row's.
read_submitter_table <- function(p) {
  x <- tryCatch(suppressWarnings(
    if (grepl("[.]xlsx?$", p, ignore.case = TRUE)) {
      as.data.frame(readxl::read_excel(p, col_names = FALSE, col_types = "text",
                                       .name_repair = "minimal"))
    } else {
      b <- readBin(p, "raw", file.size(p))
      # Windows tools often write UTF-16; a sample table is plain ASCII, so dropping the
      # zero bytes decodes it.
      if (length(b) > 1 && all(b[1:2] %in% as.raw(c(0xff, 0xfe)))) {
        b <- b[-(1:2)]
        b <- b[b != as.raw(0)]
      }
      txt <- rawToChar(b)
      # Latin-1 text read as UTF-8 would break tolower() and the join later on.
      if (!validUTF8(txt)) txt <- iconv(txt, "latin1", "UTF-8")
      txt   <- gsub("\r\n?", "\n", sub("^\ufeff", "", txt))
      lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
      n     <- function(s, l) nchar(l) - nchar(gsub(s, "", l, fixed = TRUE))
      sep   <- if (any(grepl("\t", head(lines, 10)))) "\t" else
               if (n(";", lines[1]) > n(",", lines[1])) ";" else ","
      # As many columns as the longest line, so that no line wraps onto the next.
      read.delim(text = txt, sep = sep, header = FALSE, colClasses = "character",
                 col.names = paste0("V", seq_len(1 + max(n(sep, lines)))))
    }), error = function(e) NULL)
  if (is.null(x) || nrow(x) < 2) return(NULL)
  full <- !is.na(x) & trimws(as.matrix(x)) != ""
  n    <- rowSums(head(full, 10))
  h    <- which(n >= 2 & n > median(n) / 2)[1]
  if (is.na(h)) return(NULL)
  named <- full[h, ]
  # Columns with neither a name nor a value (trailing separators) are dropped.
  keep <- named | colSums(full[-seq_len(h), , drop = FALSE]) > 0
  nm   <- ifelse(named, trimws(unlist(x[h, ])), paste0("V", seq_along(x)))[keep]
  x    <- x[-seq_len(h), keep, drop = FALSE]
  names(x)    <- make.unique(nm)
  rownames(x) <- NULL
  if (nrow(x)) x else NULL
}

# Columns of x naming the deposit's MS files: the first where >= a fifth (and >= 3) of values
# are MS-file keys, then others naming other files. None for file lists or feature tables.
msv_file_columns <- function(x, ms_keys) {
  if (ncol(x) < 2 || any(grepl(MSV_MS_EXT, names(x), ignore.case = TRUE)) ||
      any(c("#Scan#", "SpectrumID") %in% names(x)) ||
      mean(grepl("^[0-9]", names(x))) > 0.5) return(integer())
  js   <- integer()
  seen <- character()
  for (j in seq_along(x)) {
    v   <- msv_file_key(x[[j]][!is.na(x[[j]]) & nzchar(trimws(x[[j]]))])
    hit <- v %in% ms_keys
    if (!length(v) || sum(hit) < min(3, nrow(x)) || mean(hit) < 0.2) next
    if (length(js) && (anyDuplicated(v[hit]) || any(v[hit] %in% seen))) next
    js   <- c(js, j)
    seen <- c(seen, v[hit])
  }
  js
}

# Cache a MassIVE deposit's submitter sample table as <id>_metadata.tsv (path returned, or NULL);
# none found leaves <id>_metadata_all.none, a failure caches nothing. A non-empty marker holds
# the listing's warning about hidden files: delete it to search again once the package lists them.
fetch_massive_submitter_table <- function(deposit_id, cache_dir = META_DIR) {
  cached <- list.files(cache_dir, full.names = TRUE,
                       pattern = paste0("^", deposit_id, "_metadata\\.(csv|tsv)$"))
  if (length(cached)) return(cached[1])
  none <- file.path(cache_dir, paste0(deposit_id, "_metadata_all.none"))
  if (file.exists(none)) return(NULL)
  hidden <- character()
  fl <- tryCatch(withCallingHandlers(
    massive_files(deposit_id),
    warning = function(w) {
      hidden <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }), error = function(e) {
      if (is_network_error(e)) stop(e)
      NULL
    })
  if (is.null(fl)) return(NULL)                    # not listed now: try next time
  ms_keys <- msv_ms_file_keys(fl)
  cand    <- if (length(ms_keys)) msv_table_candidates(fl) else character()
  # Without readxl every Excel table would read as unreadable, and the deposit
  # would be marked as searched.
  if (any(grepl("[.]xlsx?$", cand, ignore.case = TRUE)) &&
      !requireNamespace("readxl", quietly = TRUE))
    stop("Excel sample tables need the readxl package: install.packages(\"readxl\")")
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  failed <- FALSE
  tabs <- lapply(cand, function(f) {
    download <- function() {
      unlink(list.files(tmp, full.names = TRUE))   # a partial download
      tryCatch({
        suppressWarnings(suppressMessages(MsBackendMassIVE::massive_download_file(
          deposit_id, pattern = paste0("^", gsub("([^A-Za-z0-9/_-])", "[\\1]", f), "$"),
          fileName = basename(f), path = tmp, overwrite = TRUE)))
        length(list.files(tmp)) > 0
      }, error = function(e) {
        if (is_network_error(e)) stop(e)
        FALSE
      })
    }
    # A second try, as the MassIVE FTP drops a connection now and then.
    if (!download() && !download()) {
      message("  ", deposit_id, ": could not download ", f)
      failed <<- TRUE
      return(NULL)
    }
    # The file arrives under its URL-encoded name ("new%202.txt"): take what came.
    p <- list.files(tmp, full.names = TRUE)
    x <- read_submitter_table(p[1])
    unlink(p)
    if (is.null(x)) return(NULL)
    # A table shared by several deposits (a MassiveID column): this deposit's rows.
    own <- Find(function(v) deposit_id %in% v, x)
    if (!is.null(own)) x <- x[own %in% deposit_id, , drop = FALSE]
    js <- msv_file_columns(x, ms_keys)
    if (!length(js)) return(NULL)
    message("  ", deposit_id, ": sample table ", f, " (files in '",
            paste(names(x)[js], collapse = "', '"), "')")
    # One copy of the table per file column, that column named "filename".
    dplyr::bind_rows(lapply(js, function(j) {
      names(x)[names(x) == "filename"] <- "filename_submitted"
      names(x)[j] <- "filename"
      x
    }))
  })
  # A failed download could hide the better table: keep nothing, try next time.
  if (failed) return(NULL)
  tabs <- Filter(Negate(is.null), tabs)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  if (!length(tabs)) {
    writeLines(hidden, none)
    return(NULL)
  }
  # One row per deposited file, the first table in the order above winning; rows naming files
  # never deposited would count as samples when Pan-ReDU has none.
  out <- dplyr::bind_rows(tabs)
  key <- msv_file_key(out$filename)
  out <- out[key %in% ms_keys & !duplicated(key), ]
  path <- file.path(cache_dir, paste0(deposit_id, "_metadata.tsv"))
  # Quoted, so that a stray quote or line break in a cell survives read.delim().
  write.table(out, path, sep = "\t", row.names = FALSE, quote = TRUE,
              qmethod = "double", na = "")
  path
}

# ---- extraction-set selection from the deposit inventory -----------------

# Cached names (MsBackend prefixes them) back to the bare name: the longest name of `bare_ref`
# each ends with. Here, not in ms1.R, because norm_file_key() needs it.
to_bare_name <- function(x, bare_ref) {
  if (!length(x)) return(x)
  bare_ref <- unique(bare_ref)
  # One by one, as a deposit can mix prefix forms; longest match wins, and the leading
  # "_" keeps 21_p from resolving against 1_p.
  vapply(x, function(f) {
    if (f %in% bare_ref) return(f)
    cand <- bare_ref[endsWith(f, paste0("_", bare_ref))]
    if (!length(cand)) return(f)
    cand[which.max(nchar(cand))]
  }, character(1), USE.NAMES = FALSE)
}

# Name key across sources: accession prefix and extension dropped, lower case (957.mzML -> "957").
# Pass `bare_ref` (the real bio_files) for cache-derived names: only it resolves the assay prefix.
norm_file_key <- function(x, bare_ref = NULL) {
  b <- basename(as.character(x))
  if (!is.null(bare_ref) && length(bare_ref))
    b <- to_bare_name(b, basename(as.character(bare_ref)))
  b <- sub("^(MSV[0-9]+|MTBLS[0-9]+|ST[0-9]+)[_-]", "", b)
  tolower(sub("[.][^.]*$", "", b))
}

# MS1 extraction set of a MassIVE deposit from its own file listing: files with metadata
# (meta_keys) plus hit files not declared non-biological. Returns data.frame(file, has_metadata).
select_bio_files_massive <- function(deposit_id, meta_keys = character(),
                                     hit_files = character(),
                                     nonbio_keys = character()) {
  inv <- tryCatch(massive_files(deposit_id),
                  error = function(e) {
                    if (is_network_error(e)) stop(e)
                    NULL
                  })
  if (is.null(inv)) stop("massive_list_files failed for ", deposit_id)
  paths <- if (is.data.frame(inv)) {
    fc <- intersect(c("file_name", "fileName", "name", "path", "file",
                      "spectra_file"), names(inv))[1]
    as.character(inv[[fc]])
  } else as.character(inv)

  # One copy per sample by normalised name across all folders (later updates add collections
  # elsewhere); the ccms_peak copy wins.
  ms  <- grep("[.](mzML|mzXML)$", paths, value = TRUE, ignore.case = TRUE)
  ms  <- ms[order(!grepl("^ccms_peak/", ms, ignore.case = TRUE))]
  bio <- ms[!duplicated(norm_file_key(basename(ms)))]
  with_meta <- basename(bio)[norm_file_key(basename(bio)) %in% meta_keys]
  # Plus every confirmed-hit file, whatever tree it sits in, unless the metadata declares
  # it non-biological.
  hits <- basename(as.character(hit_files))
  hits <- hits[!norm_file_key(hits) %in% nonbio_keys]
  all_files <- unique(c(with_meta, hits))

  # Pooled QCs, QC injections, blanks and standards dropped outright (EXCL_NONBIO).
  n_nonbio <- sum(is_nonbio_file(all_files))
  if (n_nonbio) {
    message("  ", deposit_id, ": dropped ", n_nonbio, " QC, blank or standard files")
    all_files <- all_files[!is_nonbio_file(all_files)]
  }

  data.frame(file         = all_files,
             has_metadata = norm_file_key(all_files) %in% meta_keys,
             stringsAsFactors = FALSE)
}

# ===== assay/technical metadata =====

# fetch_assay_metadata(dataset_ids): one row per (dataset_id, assay). Sources: MTBLS from
# mtbls_assay_data, ST from mwb_metadata()$MS_run, MSV from local TSV then Pan-ReDU fallback.

# ReDU/MIxS placeholders treated as missing.
REDU_MISSING <- c("missing value", "not applicable", "not collected",
                  "not provided", "restricted access", "na", "n/a", "",
                  "-", "unknown", "no data")

# Numeric where the text is a number, NA where it is not, decided by a pattern rather
# than by hiding as.numeric()'s warning.
as_num <- function(x) {
  x <- as.character(x)
  ok <- !is.na(x) & grepl("^ *[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)? *$", x)
  out <- rep(NA_real_, length(x))
  out[ok] <- as.numeric(x[ok])
  out
}

is_missing_val <- function(x) {
  v <- tolower(trimws(as.character(x)))
  # Pan-ReDU's UniqueSubjectID puts the accession before the recorded id, so a
  # missing one reads "MSV000094097_nan".
  v <- sub("^(msv|mtbls|st)[0-9]+_", "", v)
  is.na(x) | v %in% c(REDU_MISSING, "nan")
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
  # The assays the cached investigation file names, each read from its cached copy.
  assay_files <- grep("^a_.*\\.txt$", mtbls_inv_assay_names(ds), value = TRUE)
  if (!length(assay_files)) return(NULL)

  do.call(rbind, lapply(assay_files, function(an) {
    ad <- fetch_assay_table(ds, an)
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
  md <- mwb_metadata_cached(ds)
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
  # Submitter tables sometimes prefix the ReDU names (MSV000086158: Analysis_...,
  # ATTRIBUTE_Analysis_...); each prefixed form is matched too.
  get_col <- function(df, col) {
    if (is.null(df)) return(NULL)
    hit <- intersect(c(col, paste0("Analysis_", col), paste0("ATTRIBUTE_Analysis_", col),
                       paste0("ATTRIBUTE_", col)), names(df))
    if (!length(hit)) NULL else clean(df[[hit[1]]])
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
  rows <- lapply(unique(dataset_ids), function(ds) tryCatch({
    if      (startsWith(ds, "MSV"))   fetch_assay_metadata_massive(ds)
    else if (startsWith(ds, "MTBLS")) fetch_assay_metadata_mtbls(ds)
    else if (startsWith(ds, "ST"))    fetch_assay_metadata_mwb(ds)
    else NULL
  }, error = function(e) {
    if (is_network_error(e)) stop(e)
    message("  technical metadata not reachable for ", ds, ": ", conditionMessage(e))
    NULL
  }))
  do.call(rbind, rows)
}

# ===== biological file selection =====

META_DIR <- "dataset_metadata"
EXCL_ST  <- "(?i)^(blank|qc|pool|reference|standard|control|empty)"

# Non-biological files by name (pooled QCs, QC injections, blanks, standards), often labelled
# biological: a leading "pool" or "std<digit>", or a whole "qc" or "blank" token.
EXCL_NONBIO <- "(?i)(^pool)|((^|[_-])(start_)?qc([_-]?[0-9]+)?([_.-]|$))|((^|[_-])blank([_-]|[0-9]|$))|(^std[0-9])"
is_nonbio_file <- function(f)
  grepl(EXCL_NONBIO, sub("^(MSV[0-9]+|MTBLS[0-9]+|ST[0-9]+)[_-]", "",
                         basename(as.character(f))), perl = TRUE)

# Assay file names for a MetaboLights study from its ISA investigation file (cached locally).
mtbls_inv_assay_names <- function(ds) {
  inv_path <- file.path(META_DIR, paste0(ds, "_i_Investigation.txt"))
  if (!file.exists(inv_path)) {
    inv_url   <- paste0("https://ftp.ebi.ac.uk/pub/databases/metabolights/",
                        "studies/public/", ds, "/i_Investigation.txt")
    inv_lines <- offline_as_error(readLines(url(inv_url), warn = FALSE))
    writeLines(inv_lines, inv_path)
  } else {
    inv_lines <- readLines(inv_path, warn = FALSE)
  }
  assay_line <- inv_lines[startsWith(inv_lines, "Study Assay File Name")]
  if (!length(assay_line)) return(NULL)
  an <- trimws(unlist(strsplit(sub("^[^\t]+\t", "", assay_line[1]), "\t")))
  an[nchar(an) > 0]
}

# Fetch (or load from cache) one MetaboLights assay table; tries mtbls_assay_data() then FTP
# (FTP only when mtbls_assay_data() failed for another reason than the network).
fetch_assay_table <- function(ds, assay_name) {
  cache_path <- file.path(META_DIR, assay_name)
  if (file.exists(cache_path))
    return(read.delim(cache_path, check.names = FALSE))
  adf <- tryCatch(offline_as_error(mtbls_assay_data(ds, assay_name)),
                  error = function(e) {
                    if (is_network_error(e)) stop(e)
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
  adf <- offline_as_error(read.delim(url(assay_url), check.names = FALSE))
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
    md <- mwb_metadata_cached(ds)
    if (is.null(md) || is.null(md$MS_run) || !nrow(md$MS_run)) return(NULL)
    an_labels <- paste0(ds, "_", mwb_analysis_order(ds), ".txt")
    # A study with one analysis: every file of the deposit belongs to it.
    if (length(an_labels) == 1)
      return(stats::setNames(list(mwb_listed_files(ds)), an_labels))
    rfn_col   <- grep("RAW_FILE_NAME", names(md$sample_annotation),
                      value = TRUE)[1]
    if (is.na(rfn_col)) return(NULL)
    per_sample <- strsplit(trimws(md$sample_annotation[[rfn_col]]), "\\s+")
    # Positions only mean something when each sample lists one name per analysis.
    if (length(an_labels) > 1 && all(lengths(per_sample) <= 1)) {
      warning(ds, ": one file name per sample for ", length(an_labels),
              " analyses; its hits are left unassigned")
      return(NULL)
    }
    listed <- mwb_listed_files(ds)
    out <- list()
    for (i in seq_along(an_labels)) {
      files <- vapply(per_sample,
                      function(x) if (length(x) >= i) x[i] else NA_character_,
                      character(1))
      files <- basename(files[!is.na(files) & nzchar(files)])
      # The metadata may name a file without its extension: match on the name alone.
      stem  <- function(f) tolower(tools::file_path_sans_ext(f))
      files <- listed[match(stem(files), stem(listed))]
      files <- files[!is.na(files)]
      if (length(files)) out[[an_labels[i]]] <- unique(files)
    }
    if (length(out)) out else NULL
  } else NULL
}

# Analyses in the order a deposit lists its RAW_FILE_NAME entries, where that differs from
# MS_run (by default the i-th name is the i-th analysis).
MWB_FILE_ORDER <- list(
  ST002044 = c("AN003327", "AN003328", "AN003325", "AN003326")
)

mwb_analysis_order <- function(deposit_id) {
  ord <- MWB_FILE_ORDER[[deposit_id]]
  if (is.null(ord))
    ord <- as.character(mwb_metadata_cached(deposit_id)$MS_run$ANALYSIS_ID)
  ord
}

# File names in the deposit listing, each once; a name found in two folders cannot be told
# apart, so it is left out.
mwb_listed_files <- function(deposit_id) {
  fl <- mwb_list_files_cached(deposit_id)
  b  <- basename(vapply(fl$sample_file, utils::URLdecode, character(1),
                        USE.NAMES = FALSE))
  b[!b %in% b[duplicated(b)]]
}

# A Workbench deposit's metadata and file listing, read from the saveRDS() copy in
# dataset_metadata/ (written on first fetch); delete <id>_mwb_*.rds to fetch again.
mwb_metadata_cached <- function(deposit_id)
  mwb_cached(deposit_id, "metadata", MsBackendMetabolomicsWorkbench::mwb_metadata)

mwb_list_files_cached <- function(deposit_id)
  mwb_cached(deposit_id, "files", MsBackendMetabolomicsWorkbench::mwb_list_files)

mwb_cached <- function(deposit_id, what, fetch, tag = "mwb") {
  cache <- file.path(META_DIR, paste0(deposit_id, "_", tag, "_", what, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  x <- fetch(deposit_id)
  dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)
  # Written whole before it takes its name, so that a render stopped mid-write
  # leaves no broken copy for the next renders to read.
  saveRDS(x, paste0(cache, ".part"))
  file.rename(paste0(cache, ".part"), cache)
  x
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

# MassIVE's file listing of a deposit, cached in dataset_metadata/ like the Workbench's;
# delete <deposit>_msv_files.rds to read it anew.
massive_files <- function(deposit_id)
  mwb_cached(deposit_id, "files", MsBackendMassIVE::massive_list_files, tag = "msv")

# The deposit's .mzML/.mzXML files as listed by MassIVE: Pan-ReDU can name files never
# deposited, which make massive_sync_data_files() abort.
massive_real_ms_files <- function(deposit_id) {
  fl <- massive_files(deposit_id)
  nm <- if (is.data.frame(fl)) {
    unlist(fl[[grep("file|name|path", names(fl), ignore.case = TRUE)[1]]])
  } else fl
  bn <- basename(as.character(nm))
  bn[grepl("[.]mz(X)?ML$", bn, ignore.case = TRUE)]
}
