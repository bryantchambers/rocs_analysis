#!/usr/bin/env Rscript

source(file.path("code", "wgcna_hmm", "00_config.R"))

ensure_dirs(c(DIRS$results, DIRS$main, DIRS$extended, DIRS$reports, DIRS$logs))

scripts <- c(
  "01_data_prep.R",
  "02_wgcna_main.R",
  "03_hmm_main.R",
  "04a_tea_main.R",
  "04b_tea_variability_and_nondamaged.R",
  "08_tea_figures.R",
  "04_xrf_loco_main.R",
  "05_extended_projection.R",
  "06_report.R"
)

run_one <- function(script_name) {
  script_path <- file.path(DIRS$code, script_name)
  log_path <- file.path(DIRS$logs, paste0(tools::file_path_sans_ext(script_name), ".log"))
  log_msg("Running: Rscript ", shQuote(script_path))
  status <- system2("Rscript", args = c(script_path), stdout = log_path, stderr = log_path)
  if (!identical(status, 0L)) {
    stop(sprintf("Stage failed (%s). See log: %s", script_name, log_path))
  }
  log_msg("Completed: ", script_name, " (log: ", log_path, ")")
}

for (s in scripts) run_one(s)

write_run_metadata(
  file.path(DIRS$results, "run_workflow_metadata.tsv"),
  "run_workflow.R",
  extra = list(stages = scripts)
)

log_msg("Workflow complete.")
