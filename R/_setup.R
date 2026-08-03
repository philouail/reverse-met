# R/_setup.R — shared package loading for the reverse-metabolomics pipeline.

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
