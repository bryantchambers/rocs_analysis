#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()

progress_tsv <- Sys.getenv("BALANCED_SWEEP_PROGRESS_TSV", unset = "")
if (!nzchar(progress_tsv)) {
  progress_tsv <- file.path(BAL$logs_dir, sprintf("parameter_sweep_progress_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
progress_update <- function(step_id, status, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      lane = "balancednetwork",
      step_id = step_id,
      status = status,
      detail = detail
    ),
    progress_tsv,
    sep = "\t",
    append = file.exists(progress_tsv),
    col.names = !file.exists(progress_tsv)
  )
}

ctx <- load_balanced_context()
expr_by_core <- ctx$expr_by_core

grid <- CJ(
  power = BAL_PARAMS$grid_power,
  deepSplit = BAL_PARAMS$grid_deepSplit,
  mergeCutHeight = BAL_PARAMS$grid_mergeCutHeight,
  minModuleSize = BAL_PARAMS$grid_minModuleSize
)

run_one <- function(pars) {
  fit <- tryCatch(
    fit_consensus_from_expr(
      expr_by_core = expr_by_core,
      power = pars$power,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      minModuleSize = pars$minModuleSize
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.table(
      power = pars$power,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      minModuleSize = pars$minModuleSize,
      status = "failed",
      non_grey_modules = NA_integer_,
      grey_pct = NA_real_,
      module_size_median = NA_real_,
      error = conditionMessage(fit)
    ))
  }

  colors <- fit$colors
  tab <- table(colors)
  non_g <- names(tab)[names(tab) != "grey"]
  data.table(
    power = pars$power,
    deepSplit = pars$deepSplit,
    mergeCutHeight = pars$mergeCutHeight,
    minModuleSize = pars$minModuleSize,
    status = "ok",
    non_grey_modules = length(non_g),
    grey_pct = sum(colors == "grey") / length(colors) * 100,
    module_size_median = if (length(non_g) > 0) median(as.numeric(tab[non_g])) else NA_real_,
    error = ""
  )
}

progress_update("global", "start", sprintf("n_settings=%d", nrow(grid)))
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  step_id <- sprintf("grid_%02d", i)
  pars <- grid[i]
  if (i %% 5 == 0 || i == 1 || i == nrow(grid)) {
    log_msg(sprintf("balanced sweep %d/%d", i, nrow(grid)))
  }
  progress_update(
    step_id,
    "start",
    sprintf("power=%s deepSplit=%s mergeCutHeight=%s minModuleSize=%s", pars$power, pars$deepSplit, pars$mergeCutHeight, pars$minModuleSize)
  )
  out <- run_one(pars)
  progress_update(step_id, out$status[[1]], if (nzchar(out$error[[1]])) out$error[[1]] else "fit_complete")
  out[, error := NULL]
  out
}), fill = TRUE)

setorder(res, status, grey_pct, -non_grey_modules)
fwrite(res, file.path(BAL$qc_tables_dir, "qc_parameter_sweep_summary.tsv"), sep = "\t")

ok <- res[status == "ok"]
if (nrow(ok) > 0) {
  best <- ok[order(grey_pct, -non_grey_modules, -module_size_median)][1]
  fwrite(best, file.path(BAL$qc_tables_dir, "qc_parameter_sweep_recommended.tsv"), sep = "\t")
}

progress_update("global", "ok", sprintf("completed_n=%d", nrow(res)))
log_msg("balanced parameter sweep complete")
