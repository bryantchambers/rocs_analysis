#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()

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
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.table(
      power = pars$power,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      minModuleSize = pars$minModuleSize,
      status = "failed",
      non_grey_modules = NA_integer_,
      grey_pct = NA_real_,
      module_size_median = NA_real_
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
    module_size_median = if (length(non_g) > 0) median(as.numeric(tab[non_g])) else NA_real_
  )
}

res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  if (i %% 5 == 0 || i == 1 || i == nrow(grid)) {
    log_msg(sprintf("balanced sweep %d/%d", i, nrow(grid)))
  }
  run_one(grid[i])
}), fill = TRUE)

setorder(res, status, grey_pct, -non_grey_modules)
fwrite(res, file.path(BAL$qc_tables_dir, "qc_parameter_sweep_summary.tsv"), sep = "\t")

ok <- res[status == "ok"]
if (nrow(ok) > 0) {
  best <- ok[order(grey_pct, -non_grey_modules, -module_size_median)][1]
  fwrite(best, file.path(BAL$qc_tables_dir, "qc_parameter_sweep_recommended.tsv"), sep = "\t")
}

log_msg("balanced parameter sweep complete")
