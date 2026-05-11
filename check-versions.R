# check-versions.R
# Source this at the top of any script/document that uses RuSirius.
# Stops with install instructions if a required version is not met.

required <- list(
  RuSirius         = list(version = "1.0.1",
                          install = 'remotes::install_github("RforMassSpectrometry/RuSirius@v1.0.1")'),
  MsBackendMassIVE = list(version = "0.99.0",
                          install = 'remotes::install_github("RforMassSpectrometry/MsBackendMassIVE@v0.99.0")')
)

for (pkg in names(required)) {
  req  <- required[[pkg]]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is not installed.\nRun:\n  ", req$install, call. = FALSE)
  }
  got <- as.character(packageVersion(pkg))
  if (got != req$version) {
    stop(
      pkg, " version ", got, " is installed, but ", req$version, " is required.\n",
      "Run:\n  ", req$install,
      call. = FALSE
    )
  }
}

# SIRIUS binary — only major.minor must match (6.4.x)
SIRIUS_REQUIRED <- "6.4"

sirius_bin <- Sys.which("sirius")
if (nchar(sirius_bin) == 0) {
  stop(
    "SIRIUS binary not found on PATH.\n",
    "Download SIRIUS 6.4 from https://github.com/boecker-lab/sirius/releases\n",
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
    "Expected SIRIUS ", SIRIUS_REQUIRED, ".x.",
    call. = FALSE
  )
}
sirius_ver      <- ver_match[1]
sirius_maj_min  <- paste(strsplit(sirius_ver, "\\.")[[1]][1:2], collapse = ".")
if (sirius_maj_min != SIRIUS_REQUIRED) {
  stop(
    "SIRIUS ", sirius_ver, " found, but ", SIRIUS_REQUIRED, ".x is required.\n",
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
