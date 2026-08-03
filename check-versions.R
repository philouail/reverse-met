# check-versions.R — source at the top of any pipeline script/document.
# Stops with install instructions if a required version is not met.
# `type`: "exact" (must match — RuSirius is JOSS-pinned) or "min" (>= version).

required <- list(
  RuSirius                     = list(version = "1.0.4", type = "exact",
                                      install = 'remotes::install_github("RforMassSpectrometry/RuSirius@v1.0.4")'),
  # gabri branches: bake the assay index into the BiocFileCache filename so a
  # deposit's pos/neg files with identical basenames no longer collide (loading
  # the wrong polarity). Required for correct per-assay downloads.
  MsBackendMassIVE             = list(version = "0.99.1", type = "min",
                                      install = 'remotes::install_github("RforMassSpectrometry/MsBackendMassIVE@gabri")'),
  MsBackendMetaboLights        = list(version = "1.7.4",  type = "min",
                                      install = 'remotes::install_github("RforMassSpectrometry/MsBackendMetaboLights@gabri")'),
  MsBackendMetabolomicsWorkbench = list(version = "0.1.3", type = "min",
                                      install = 'BiocManager::install("MsBackendMetabolomicsWorkbench")')
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

# SIRIUS binary — minimum version 6.3 (6.3.x and 6.4.x both supported)
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
  "Version check passed — ",
  paste(
    paste0(names(required), " ", sapply(required, `[[`, "version")),
    collapse = ", "
  ),
  ", SIRIUS ", sirius_ver
)
