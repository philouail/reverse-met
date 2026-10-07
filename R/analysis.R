# R/analysis.R: helpers shared by the two step-6 notebooks (06-analysis.qmd).
# Source after R/metadata.R, which holds norm_file_key().

# Does a metadata source cover each file? Step 5's per-assay flag, by assay and file name;
# stops when an assay has no file list (step 5 has not finished with it).
has_metadata <- function(pf, bio_meta) {
  miss <- setdiff(unique(pf$dataset_id), names(bio_meta))
  if (length(miss)) stop("no file list from step 5 for ", paste(miss, collapse = ", "))
  key <- unlist(lapply(names(bio_meta), function(ds) {
    v <- bio_meta[[ds]]
    stats::setNames(as.logical(v), paste0(ds, "|", basename(names(v))))
  }))
  h <- key[paste0(pf$dataset_id, "|", basename(as.character(pf$file)))]
  unname(!is.na(h) & h)
}

# TRUE where the file gave a SIRIUS-confirmed spectrum of each target (step 5's MS2 rule).
# Confirmed spectra use cache file names, so keys are resolved against the real file list.
ms2_confirmed <- function(pf, conf, bio_files) {
  miss <- setdiff(unique(pf$dataset_id), names(bio_files))
  if (length(miss)) stop("no file list from step 5 for ", paste(miss, collapse = ", "))
  out <- logical(nrow(pf))
  for (d in unique(pf$dataset_id)) {
    i  <- pf$dataset_id == d
    fb <- basename(as.character(bio_files[[d]]))
    k  <- norm_file_key(pf$file[i])
    has <- function(role) k %in% norm_file_key(
      conf$source_file[conf$dataset_id == d & conf$role == role], fb)
    out[i] <- has("parent") & has("metabolite")
  }
  out
}

# Groups of assays measuring the same samples (shared sample names within a deposit, or a
# panel declared in related_deposits.csv), each named after its first assay id.
same_sample_groups <- function(pf, datasets, related) {
  ids  <- sort(unique(pf$dataset_id), method = "radix")
  miss <- setdiff(ids, datasets$dataset_id)
  if (length(miss)) stop("not among the confirmed datasets: ", paste(miss, collapse = ", "))
  dep <- datasets$deposit_id[match(ids, datasets$dataset_id)]
  sn  <- lapply(ids, function(d) {
    s <- as.character(pf$sample_name[pf$dataset_id == d])
    unique(s[!is.na(s) & nzchar(s)])
  })

  # Pairs of assays measuring the same samples: within one deposit, at least
  # half the sample names of the smaller one are shared.
  pairs <- list()
  for (i in seq_along(ids)) for (j in seq_along(ids)) {
    if (i >= j || dep[i] != dep[j]) next
    shared <- length(intersect(sn[[i]], sn[[j]]))
    if (shared > 0 && shared >= 0.5 * min(length(sn[[i]]), length(sn[[j]])))
      pairs[[length(pairs) + 1]] <- c(ids[i], ids[j])
  }
  # Panels are declared by accession; one split into several assays is ambiguous, so it
  # stops. A panel with a side outside the corpus is skipped.
  to_dataset <- function(x) {
    d <- ids[ids %in% x | dep %in% x]
    if (length(d) > 1) stop(x, " has several assays: declare its panel by dataset id")
    d
  }
  panel <- related[related$relation == "same_samples_two_assays", , drop = FALSE]
  for (r in seq_len(nrow(panel))) {
    p <- c(to_dataset(panel$dataset_id[r]), to_dataset(panel$related_to[r]))
    if (length(p) == 2) pairs[[length(pairs) + 1]] <- p
  }

  # Join the pairs into groups.
  grp <- stats::setNames(ids, ids)
  for (p in pairs) grp[grp == grp[[p[2]]]] <- grp[[p[1]]]
  first_id <- vapply(grp, function(g) sort(ids[grp == g], method = "radix")[1], "")
  data.frame(dataset_id = ids, group = unname(first_id), stringsAsFactors = FALSE)
}

# Keeps one assay per same-sample group: most files with both targets by MS1, then by MS2,
# then the first id. One row per assay: dataset_id, group, files, ms1_both, ms2_both, kept.
keep_one_assay <- function(pf, datasets, related, conf, bio_files) {
  g <- same_sample_groups(pf, datasets, related)
  f <- factor(pf$dataset_id, levels = g$dataset_id)
  g$files    <- as.integer(table(f))
  g$ms1_both <- as.integer(tapply(pf$par_detected %in% TRUE & pf$met_detected %in% TRUE, f, sum))
  g$ms2_both <- as.integer(tapply(ms2_confirmed(pf, conf, bio_files), f, sum))
  g <- g[order(g$group, -g$ms1_both, -g$ms2_both, g$dataset_id, method = "radix"), ]
  g$kept <- !duplicated(g$group)
  rownames(g) <- NULL
  g
}

# Empty text as NA.
nz <- function(v) ifelse(is.na(v) | !nzchar(as.character(v)), NA_character_, as.character(v))

# Sex is deposited as F/M by some deposits and female/male by the rest.
harmonise_sex <- function(x) {
  s <- tolower(trimws(as.character(x)))
  dplyr::case_when(startsWith(s, "f") ~ "female",
                   startsWith(s, "m") ~ "male",
                   TRUE ~ NA_character_)
}

# Deposits name the same matrix differently (blood_plasma and "blood plasma",
# UANoAdditive_YellowTop and urine, fecal and feces): one name per matrix.
canon_site <- function(x) {
  y <- tolower(trimws(gsub("[_.]+", " ", as.character(x))))
  dplyr::case_when(
    is.na(y) | !nzchar(y)               ~ NA_character_,
    grepl("plasma", y)                  ~ "plasma",
    grepl("serum", y)                   ~ "serum",
    grepl("^ua|urine", y)               ~ "urine",
    grepl("fecal|feces|faec|stool", y)  ~ "faeces",
    grepl("forehead|arm skin|^skin", y) ~ "skin",
    grepl("saliva", y)                  ~ "saliva",
    grepl("oral cavity|mouth", y)       ~ "oral cavity",
    grepl("milk", y)                    ~ "milk",
    grepl("nasal", y)                   ~ "nasal cavity",
    grepl("vagina", y)                  ~ "vagina",
    grepl("anal", y)                    ~ "anal region",
    TRUE                                ~ y)
}

# One unit per subject and matrix (per file where no subject is named): l2 is the median ratio,
# carried columns their one value or NA; disagreeing units/fields go in attr "disagree".
unitise <- function(x, fields = character(0)) {
  x$unit <- ifelse(is.na(x$subj), paste("file", x$study, x$fk),
                   paste(x$study, x$subj, x$site))
  carry <- union(intersect(c("sex2", "agen", "health"), names(x)), fields)
  if ("health" %in% names(x)) x$health[x$health %in% "not recorded"] <- NA
  one_value <- function(v) {
    u <- unique(v[!is.na(v)])
    if (length(u) == 1) u else v[NA_integer_]
  }
  out <- x |>
    dplyr::group_by(unit) |>
    dplyr::summarise(
      study  = dplyr::first(study),
      kind   = if (is.na(dplyr::first(subj))) "file" else "subject",
      site   = dplyr::first(site),
      files  = dplyr::n(),
      either = any(par_detected %in% TRUE | met_detected %in% TRUE),
      l2     = if (any(is.finite(l2))) stats::median(l2[is.finite(l2)]) else NA_real_,
      dplyr::across(dplyr::all_of(carry), one_value),
      .groups = "drop") |>
    as.data.frame()
  if ("health" %in% names(out)) out$health[is.na(out$health)] <- "not recorded"
  # The units whose files record more than one value of a column.
  dis <- data.frame(unit = character(0), field = character(0))
  for (f in carry) {
    n <- tapply(x[[f]], x$unit, function(v) length(unique(v[!is.na(v)])))
    if (any(n > 1)) dis <- rbind(dis, data.frame(unit = names(n)[n > 1], field = f))
  }
  attr(out, "disagree") <- dis
  out
}

# A unitise() result's disagreements as text, "unit (fields); ...", fields named by `labels`;
# "" where none.
disagree_text <- function(u, labels = character(0)) {
  d <- attr(u, "disagree")
  if (nrow(d) == 0) return("")
  field <- ifelse(d$field %in% names(labels), labels[d$field], d$field)
  f <- tapply(field, d$unit, paste, collapse = ", ")
  paste0(names(f), " (", f, ")", collapse = "; ")
}
