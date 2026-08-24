# R/masst.R — read and pre-process StructureMASST exports

ADDUCT_CANDIDATES <- c(
  "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+",
  "[M-H]-", "[M+Cl]-"
)

# Read one MASST CSV, label its role, and split the USI into file_path/scan_str.
read_masst <- function(path, role) {
  x <- read.csv(path, check.names = FALSE)
  for (j in seq_along(x))
    if (is.character(x[[j]])) x[[j]][is_missing_val(x[[j]])] <- NA
  parts       <- strsplit(x$USI, ":", fixed = TRUE)
  x$file_path <- vapply(parts, `[`, character(1), 3)
  x$scan_str  <- vapply(parts, \(p) p[length(p)], character(1))
  x$role      <- role
  x
}

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

# ===== compound pair (was compounds.R) =====

# R/compounds.R — single source of truth for the parent/metabolite pair under
# study. To swap in a different probe pair, edit this file (and re-run StructureMASST).

compound_spec <- data.frame(
  role    = c("parent",       "metabolite"),
  name    = c("omeprazole",   "5-hydroxyomeprazole"),
  formula = c("C17H19N3O3S",  "C17H19N3O4S"),
  mz_pos  = c(346.1220,       362.1169),
  inchikey_block1 = c("SUBDBMMJDZJVOS", "CMZHQFXXAAIBKE"),
  stringsAsFactors = FALSE
)

# CYP2C19 activity order, low -> high. The ENTIRE validation is a regression on this
# ordering, so its direction sets the sign of the headline result. It was defined twice,
# verbatim, in 03 and 04; reordering one copy would have flipped a slope with nothing
# raising an error. One definition, here with the other domain constants.
PHENO <- c("Poor Metabolizer", "Intermediate Metabolizer", "Normal Metabolizer",
           "Rapid Metabolizer", "Ultrarapid Metabolizer")

# role -> formula, the form most call sites want.
role_formula <- setNames(compound_spec$formula, compound_spec$role)
