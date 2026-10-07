# R/setup.R: sourced by every notebook. Checks pinned versions and the SIRIUS binary,
# loads the shared packages and defines the helpers that name the targets.

# ===== versions =====
# Stops with install instructions if a pinned package or the SIRIUS binary is missing
# or off-version. `type`: "exact" (must match; RuSirius is JOSS-pinned) or "min" (>= version).
check_versions <- function() {
  required <- list(
    RuSirius                     = list(version = "1.0.4", type = "exact",
                                        install = 'remotes::install_github("RforMassSpectrometry/RuSirius@v1.0.4")'),
    # gabri branches put the assay index in the cache filename, so a deposit's pos/neg
    # files with the same basename do not collide (loading the wrong polarity).
    MsBackendMassIVE             = list(version = "0.99.1", type = "min",
                                        install = 'remotes::install_github("RforMassSpectrometry/MsBackendMassIVE@gabri")'),
    MsBackendMetaboLights        = list(version = "1.7.4",  type = "min",
                                        install = 'remotes::install_github("RforMassSpectrometry/MsBackendMetaboLights@gabri")'),
    MsBackendMetabolomicsWorkbench = list(version = "0.99.0", type = "min",
                                        install = 'remotes::install_github("RforMassSpectrometry/MsBackendMetabolomicsWorkbench")'),
    # The MS1 extraction runs on Chromatograms; the results were produced with 1.3.3.
    Chromatograms                = list(version = "1.3.3",  type = "min",
                                        install = 'BiocManager::install("Chromatograms", update = FALSE, ask = FALSE)')
  )

  for (pkg in names(required)) {
    req <- required[[pkg]]
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg, " is not installed.\nRun:\n  ", req$install, call. = FALSE)
    }
    got <- as.character(packageVersion(pkg))
    ok  <- if (req$type == "exact") got == req$version
           else                     utils::compareVersion(got, req$version) >= 0
    if (!ok) {
      op <- if (req$type == "exact") "is required" else "or newer is required"
      stop(
        pkg, " version ", got, " is installed, but ", req$version, " ", op, ".\n",
        "Run:\n  ", req$install,
        call. = FALSE
      )
    }
  }

  # SIRIUS binary: minimum version 6.3 (6.3.x and 6.4.x both supported)
  SIRIUS_MIN <- "6.3"

  sirius_bin <- Sys.which("sirius")
  if (nchar(sirius_bin) == 0) {
    stop(
      "SIRIUS binary not found on PATH.\n",
      "Download SIRIUS >= ", SIRIUS_MIN, " from ",
      "https://github.com/boecker-lab/sirius/releases\n",
      "and make sure it is on your PATH.",
      call. = FALSE
    )
  }

  raw <- tryCatch(
    system2(sirius_bin, "--version", stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  ver_match <- regmatches(raw, regexpr("[0-9]+\\.[0-9]+(?:\\.[0-9]+)?", raw))
  if (!length(ver_match)) {
    stop(
      "Could not parse SIRIUS version from '", sirius_bin, " --version'.\n",
      "Expected SIRIUS >= ", SIRIUS_MIN, ".",
      call. = FALSE
    )
  }
  sirius_ver <- ver_match[1]
  if (utils::compareVersion(sirius_ver, SIRIUS_MIN) < 0) {
    stop(
      "SIRIUS ", sirius_ver, " found, but >= ", SIRIUS_MIN, " is required.\n",
      "Download from https://github.com/boecker-lab/sirius/releases",
      call. = FALSE
    )
  }

  message(
    "Version check passed: ",
    paste(
      paste0(names(required), " ", sapply(required, `[[`, "version")),
      collapse = ", "
    ),
    ", SIRIUS ", sirius_ver
  )
}
check_versions()


suppressPackageStartupMessages({
  library(Spectra)
  library(MsBackendMassIVE)
  library(MsBackendMetaboLights)
  library(MsBackendMetabolomicsWorkbench)
  library(xcms)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MetaboCoreUtils)
  library(Chromatograms)
  library(MsQuality)
  library(MsCoreUtils)
  library(RuSirius)
  library(jsonlite)   # Pan-ReDU JSON parsing (R/metadata.R)
  library(curl)       # Pan-ReDU HTTPS fetch with SSL bypass (R/metadata.R)
})

dir.create("artifacts", showWarnings = FALSE)

# ===== the targets =====
# compound.csv row order is the target order. Per-target columns carry a prefix: par
# (parent), met (metabolite), else the role itself (par_auc, mz_met, sulfone_auc).

target_prefix <- function(role) {
  p <- as.character(role)
  p[p %in% "parent"]     <- "par"
  p[p %in% "metabolite"] <- "met"
  p
}

role_of_prefix <- function(prefix) {
  r <- as.character(prefix)
  r[r %in% "par"] <- "parent"
  r[r %in% "met"] <- "metabolite"
  r
}

# How a step summary names the targets: "the target", "both targets",
# "all 3 targets".
targets_phrase <- function(n) {
  if (n == 1) "the target" else if (n == 2) "both targets" else paste("all", n, "targets")
}

# ... and an assay that confirms some of them but not all: "one target" when
# there are two.
some_targets_phrase <- function(n) {
  if (n == 2) "one target" else "some but not all targets"
}
