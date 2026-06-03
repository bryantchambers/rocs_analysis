#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

OUT_BASE <- here("networkQC", "results")
OUT_TABLE <- file.path(OUT_BASE, "tables")
OUT_FIG <- file.path(OUT_BASE, "figures")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

req <- c(
  file.path(RESULTS$stage1, "prokaryotes_vst.rds"),
  file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"),
  file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"),
  file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv")
)
missing <- req[!file.exists(req)]
if (length(missing) > 0) {
  stop("Missing required inputs:\n", paste0(" - ", missing, collapse = "\n"))
}

params <- data.table(
  param = c(
    "seed", "stage1_cores", "validation_core",
    "sweep_powers", "sweep_deepSplit", "sweep_mergeCutHeight", "sweep_minModuleSize",
    "age_grid_points", "leiden_resolution", "leiden_tom_threshold"
  ),
  value = c(
    PARAMS$seed,
    paste(PARAMS$stage1_cores, collapse = ","),
    PARAMS$validation_core,
    "10,12,14,16,18,20",
    "1,2,3",
    "0.10,0.15,0.20,0.25",
    "20,30,40",
    "100",
    "0.05,0.1,0.2,0.4",
    "0.05"
  )
)
fwrite(params, file.path(OUT_TABLE, "network_qc_parameters.tsv"), sep = "\t")
message("[networkQC] config check passed")
