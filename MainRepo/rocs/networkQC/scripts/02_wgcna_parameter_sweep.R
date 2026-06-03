#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

OUT_TABLE <- here("networkQC", "results", "tables")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

multiExpr <- lapply(PARAMS$stage1_cores, function(core_id) list(data = expr_by_core[[core_id]]))
names(multiExpr) <- PARAMS$stage1_cores

grid <- CJ(
  power = c(12, 16, 20),
  deepSplit = c(1, 2, 3),
  mergeCutHeight = c(0.10, 0.15, 0.20, 0.25),
  minModuleSize = c(20, 30)
)

run_one <- function(pars) {
  fit <- tryCatch(
    blockwiseConsensusModules(
      multiExpr = multiExpr,
      power = pars$power,
      networkType = "signed",
      corType = "pearson",
      maxBlockSize = 5000,
      minModuleSize = pars$minModuleSize,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      numericLabels = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.table(
      power = pars$power, deepSplit = pars$deepSplit, mergeCutHeight = pars$mergeCutHeight, minModuleSize = pars$minModuleSize,
      status = "failed", non_grey_modules = NA_integer_, grey_pct = NA_real_, module_size_median = NA_real_
    ))
  }
  colors <- fit$colors
  tab <- table(colors)
  non_g <- names(tab)[names(tab) != "grey"]
  data.table(
    power = pars$power, deepSplit = pars$deepSplit, mergeCutHeight = pars$mergeCutHeight, minModuleSize = pars$minModuleSize,
    status = "ok",
    non_grey_modules = length(non_g),
    grey_pct = sum(colors == "grey") / length(colors) * 100,
    module_size_median = if (length(non_g) > 0) median(as.numeric(tab[non_g])) else NA_real_
  )
}

res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  if (i %% 5 == 0 || i == 1 || i == nrow(grid)) {
    log_msg(sprintf("Sweep %d/%d", i, nrow(grid)))
  }
  run_one(grid[i])
}), fill = TRUE)

setorder(res, status, grey_pct, -non_grey_modules)
fwrite(res, file.path(OUT_TABLE, "qc_parameter_sweep_summary.tsv"), sep = "\t")

ok <- res[status == "ok"]
if (nrow(ok) > 0) {
  best <- ok[order(grey_pct, -non_grey_modules, -module_size_median)][1]
  fwrite(best, file.path(OUT_TABLE, "qc_parameter_sweep_recommended.tsv"), sep = "\t")
}

log_msg("Parameter sweep complete")

