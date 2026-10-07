# R/masst.R: adduct inference and the target constants

ADDUCT_CANDIDATES <- c(
  "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+",
  "[M-H]-", "[M+Cl]-"
)

# Infer the adduct: match observed precursor m/z to theoretical adduct m/z,
# returning the closest within ppm tolerance (else NA).
infer_adduct <- function(obs_mz, formula, ppm = 25,
                         candidates = ADDUCT_CANDIDATES) {
  if (is.na(obs_mz) || is.na(formula)) return(NA_character_)
  mass     <- calculateMass(formula)
  theo     <- vapply(candidates,
                     function(a) mass2mz(mass, adduct = a)[[1]],
                     numeric(1))
  diff_ppm <- abs(obs_mz - theo) / theo * 1e6
  best     <- which.min(diff_ppm)
  if (diff_ppm[best] > ppm) NA_character_ else candidates[best]
}

# ===== the targets =====

# The targets come from application/<name>/compound.csv (COMPOUND_SPEC overrides; default
# omeprazole). Optional isobar/isobar_at name a same-m/z compound and the target it co-elutes with.

# Roles a target may not take: the par/met prefixes, names already used by fields or
# labels, and Windows device names (a role names a directory).
RESERVED_ROLES <- c("par", "met", "all", "none", "both", "neither", "metab", "total", "ppm",
                    "con", "prn", "aux", "nul", paste0("com", 1:9), paste0("lpt", 1:9))

load_compound_spec <- function(path = Sys.getenv("COMPOUND_SPEC", "compound.csv")) {
  if (!nzchar(path) || !file.exists(path))
    return(data.frame(
      role    = c("parent",       "metabolite"),
      name    = c("omeprazole",   "5-hydroxyomeprazole"),
      formula = c("C17H19N3O3S",  "C17H19N3O4S"),
      mz_pos  = c(346.1220,       362.1169),
      inchikey_block1 = c("SUBDBMMJDZJVOS", "CMZHQFXXAAIBKE"),
      isobar  = c("hydroxyomeprazole sulfide", "omeprazole sulfone"),
      stringsAsFactors = FALSE))
  s <- read.csv(path, stringsAsFactors = FALSE)
  need <- c("role", "name", "formula", "mz_pos", "inchikey_block1")
  miss <- setdiff(need, names(s))
  if (length(miss))
    stop(path, ": missing column(s) ", paste(miss, collapse = ", "))
  if (!nrow(s))
    stop(path, ": no target; give one row per target")
  bad <- s$role[!grepl("^[a-z][a-z0-9]*$", s$role) | s$role %in% RESERVED_ROLES]
  if (length(bad))
    stop(path, ": role ", paste(bad, collapse = ", "), " is not allowed; a role is ",
         "a lowercase letter followed by lowercase letters or digits, and none of ",
         paste(RESERVED_ROLES, collapse = ", "))
  if (anyDuplicated(s$role))
    stop(path, ": role ", s$role[duplicated(s$role)][1], " is given twice")
  if (anyDuplicated(s$inchikey_block1)) {
    k <- s$inchikey_block1[duplicated(s$inchikey_block1)][1]
    stop(path, ": targets ", paste(s$role[s$inchikey_block1 == k], collapse = " and "),
         " carry the same InChIKey block ", k)
  }
  blank <- function(x) ifelse(is.na(x), "", trimws(as.character(x)))
  s$isobar <- if ("isobar" %in% names(s)) blank(s$isobar) else ""
  # isobar_at is kept only where compound.csv has it.
  has_at <- "isobar_at" %in% names(s)
  if (has_at) {
    s$isobar_at <- blank(s$isobar_at)
    lone <- nzchar(s$isobar_at) & !nzchar(s$isobar)
    if (any(lone))
      stop(path, ": target ", s$role[lone][1], " has an isobar_at but no isobar")
  }
  ref <- isobar_reference(s)
  for (i in which(nzchar(s$isobar))) {
    if (!nzchar(ref[i]))
      stop(path, ": target ", s$role[i], " declares an isobar but no isobar_at; name ",
           "the target whose retention time the isobar shares in an isobar_at column ",
           "(only with two targets is it the other one by default)")
    if (!ref[i] %in% setdiff(s$role, s$role[i]))
      stop(path, ": target ", s$role[i], " has isobar_at '", ref[i],
           "', which is not another target's role")
  }
  s[, c(need, "isobar", if (has_at) "isobar_at")]
}

# Per target, the role whose retention time its isobar shares (isobar_at, else the other of
# two targets); NA where no isobar is declared, "" where one is but no target is named.
isobar_reference <- function(spec = compound_spec) {
  n   <- nrow(spec)
  iso <- if ("isobar" %in% names(spec)) spec$isobar else rep("", n)
  at  <- if ("isobar_at" %in% names(spec)) spec$isobar_at else rep("", n)
  if (n == 2) at[!nzchar(at)] <- rev(spec$role)[!nzchar(at)]
  ifelse(!is.na(iso) & nzchar(iso), at, NA_character_)
}

compound_spec <- load_compound_spec()

# CYP2C19 activity order, low -> high. The validation regresses on this order, so its
# direction sets the sign of the headline result.
PHENO <- c("Poor Metabolizer", "Intermediate Metabolizer", "Normal Metabolizer",
           "Rapid Metabolizer", "Ultrarapid Metabolizer")

# role -> formula, the form most call sites want.
role_formula <- setNames(compound_spec$formula, compound_spec$role)

# Where a dataset's SIRIUS results are: the Workbench notebook runs per dataset,
# the MassIVE and MetaboLights notebooks per deposit.
results_dir_for <- function(dataset_id, deposit_id, results_dir = "results") {
  d <- file.path(results_dir, dataset_id)
  if (dir.exists(d)) d else file.path(results_dir, deposit_id)
}

# ===== what the structure ranking prefers =====
# One row per confirmed spectrum: the target's rank among its own formula's candidates, its
# within-spectrum score lead (if first) or deficit, and the winner (top_name, top_key).
structure_ranks <- function(conf, spec = compound_spec, results_dir = "results") {
  if (!dir.exists(results_dir)) return(NULL)
  per_role <- lapply(spec$role, function(role) {
    key  <- spec$inchikey_block1[spec$role == role]
    form <- spec$formula[spec$role == role]
    cr   <- conf[conf$role %in% role, , drop = FALSE]
    per_dep <- lapply(unique(as.character(cr$dataset_id)), function(ds) {
      dep <- as.character(cr$deposit_id[cr$dataset_id == ds][1])
      f <- file.path(results_dir_for(ds, dep, results_dir), paste0(role, "_results.csv"))
      if (!file.exists(f)) return(NULL)
      d <- read.csv(f, stringsAsFactors = FALSE)
      keep <- d$xcms_fts %in% cr$xcms_fts[as.character(cr$dataset_id) == ds] &
              d$molecularFormula %in% form & is.finite(d$csiScore)
      d <- d[keep, , drop = FALSE]
      if (!nrow(d)) return(NULL)
      d$.tgt <- startsWith(as.character(d$inchiKey), key)
      rows <- lapply(split(d, d$xcms_fts), function(x) {
        x <- x[order(-x$csiScore), , drop = FALSE]
        w <- which(x$.tgt)
        if (!length(w)) return(NULL)
        nm <- as.character(x$structureName[1])
        data.frame(
          role = role, deposit_id = dep, feature = as.character(x$xcms_fts[1]),
          candidates = nrow(x), rank = w[1],
          top_name = if (is.na(nm) || !nzchar(nm)) substr(as.character(x$inchiKey[1]), 1, 14) else nm,
          top_key  = substr(as.character(x$inchiKey[1]), 1, 14),
          lead     = if (w[1] == 1L && nrow(x) > 1L) x$csiScore[1] - x$csiScore[2] else NA_real_,
          deficit  = if (w[1] == 1L) NA_real_ else x$csiScore[1] - x$csiScore[w[1]],
          stringsAsFactors = FALSE)
      })
      do.call(rbind, rows)
    })
    do.call(rbind, per_dep)
  })
  out <- do.call(rbind, per_role)
  if (is.null(out) || !nrow(out)) return(NULL)
  rownames(out) <- NULL
  out
}
