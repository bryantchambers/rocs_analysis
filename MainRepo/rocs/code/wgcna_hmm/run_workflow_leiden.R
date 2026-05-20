#!/usr/bin/env Rscript

source(file.path("code", "wgcna_hmm", "00_config.R"))

leiden_results_root <- env_string("WGCNA_HMM_RESULTS_ROOT", default = file.path(REPO_ROOT, "results", "wgcna_hmm_leiden"))
leiden_resolution <- env_numeric("WGCNA_HMM_LEIDEN_RESOLUTION", default = 1.0)
leiden_resolution_strict <- env_flag("WGCNA_HMM_LEIDEN_RESOLUTION_STRICT", default = FALSE)

leiden_env <- c(
  paste0("WGCNA_HMM_RESULTS_ROOT=", leiden_results_root),
  "WGCNA_HMM_MODULE_METHOD=leiden",
  paste0("WGCNA_HMM_LEIDEN_RESOLUTION=", leiden_resolution),
  paste0("WGCNA_HMM_LEIDEN_RESOLUTION_STRICT=", tolower(as.character(leiden_resolution_strict)))
)

ensure_dirs(c(leiden_results_root, file.path(leiden_results_root, "logs")))

scripts <- c(
  "01_data_prep.R",
  "02_wgcna_main.R",
  "03_hmm_main.R",
  "04a_tea_main.R",
  "04b_tea_variability_and_nondamaged.R",
  "08_tea_figures.R",
  "04_xrf_loco_main.R",
  "05_extended_projection.R",
  "06_report.R",
  "09_compare_baseline_vs_leiden.R"
)

run_one <- function(script_name) {
  script_path <- file.path(DIRS$code, script_name)
  base_name <- tools::file_path_sans_ext(script_name)
  log_path <- file.path(leiden_results_root, "logs", paste0(base_name, ".log"))
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

  log_msg("[LEIDEN] Running: Rscript ", shQuote(script_path))
  status <- system2("Rscript", args = c(script_path), env = leiden_env, stdout = log_path, stderr = log_path)
  if (!identical(status, 0L)) {
    stop(sprintf("Leiden stage failed (%s). See log: %s", script_name, log_path))
  }
  log_msg("[LEIDEN] Completed: ", script_name, " (log: ", log_path, ")")
}

for (s in scripts) run_one(s)

write_run_metadata(
  file.path(leiden_results_root, "run_workflow_leiden_metadata.tsv"),
  "run_workflow_leiden.R",
  extra = list(
    stages = scripts,
    leiden_resolution_requested = leiden_resolution,
    leiden_resolution_strict = leiden_resolution_strict,
    results_root = leiden_results_root
  )
)

log_msg("Leiden alternative workflow complete. Results: ", leiden_results_root)
