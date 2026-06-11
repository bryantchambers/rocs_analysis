#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mode <- Sys.getenv("WGCNA_HMM_WGCNA_MODE", unset = "build")
if (length(args) > 0) {
  mode_arg <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
  if (!is.na(mode_arg) && nzchar(mode_arg)) mode <- mode_arg
}
if (!mode %in% c("build", "final")) {
  stop("Invalid --mode value: ", mode, ". Use --mode=build or --mode=final.")
}
Sys.setenv(WGCNA_HMM_WGCNA_MODE = mode)

source(file.path("code", "wgcna_hmm", "00_config.R"))

ensure_dirs(c(DIRS$results, DIRS$main, DIRS$extended, DIRS$reports, DIRS$qc, DIRS$wgcna_stability, DIRS$logs))

progress_tsv <- file.path(DIRS$logs, "workflow_progress.tsv")
if (file.exists(progress_tsv)) file.remove(progress_tsv)
progress_update <- function(stage, status, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stage = stage,
      status = status,
      detail = detail
    ),
    progress_tsv,
    sep = "\t",
    append = file.exists(progress_tsv),
    col.names = !file.exists(progress_tsv)
  )
}

scripts <- c(
  "01_data_prep.R",
  "02_wgcna_main.R",
  "02b_wgcna_stability.R",
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
  progress_update(script_name, "start", paste0("mode=", mode))
  script_args <- c(script_path)
  if (script_name %in% c("02_wgcna_main.R", "02b_wgcna_stability.R")) {
    script_args <- c(script_args, paste0("--mode=", mode))
  }
  status <- system2("Rscript", args = script_args, stdout = log_path, stderr = log_path)
  if (!identical(status, 0L)) {
    progress_update(script_name, "failed", log_path)
    stop(sprintf("Stage failed (%s). See log: %s", script_name, log_path))
  }
  progress_update(script_name, "ok", log_path)
  log_msg("Completed: ", script_name, " (log: ", log_path, ")")
}

progress_update("workflow", "start", paste(scripts, collapse = ","))
for (s in scripts) run_one(s)
progress_update("workflow", "ok", paste0("mode=", mode))

write_run_metadata(
  file.path(DIRS$results, "run_workflow_metadata.tsv"),
  "run_workflow.R",
  extra = list(
    stages = scripts,
    wgcna_mode = mode,
    wgcna_input_strategy = PARAMS$wgcna_input_strategy,
    wgcna_profile = PARAMS$wgcna_profile,
    hmm_input_balancing_status = PARAMS$hmm_input_balancing_status
  )
)

log_msg("Workflow complete.")
