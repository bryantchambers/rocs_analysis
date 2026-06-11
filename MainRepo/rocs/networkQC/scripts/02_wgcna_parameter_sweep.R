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
OUT_LOG <- here("networkQC", "results", "logs")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOG, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
progress_tsv <- Sys.getenv("NETWORKQC_SWEEP_PROGRESS_TSV", unset = "")
if (!nzchar(progress_tsv)) {
  progress_tsv <- file.path(OUT_LOG, sprintf("parameter_sweep_progress_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
progress_update <- function(step_id, status, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      lane = "networkQC",
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
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(data.table(
      power = pars$power, deepSplit = pars$deepSplit, mergeCutHeight = pars$mergeCutHeight, minModuleSize = pars$minModuleSize,
      status = "failed", non_grey_modules = NA_integer_, grey_pct = NA_real_, module_size_median = NA_real_,
      error = conditionMessage(fit)
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
    module_size_median = if (length(non_g) > 0) median(as.numeric(tab[non_g])) else NA_real_,
    error = ""
  )
}

progress_update("global", "start", sprintf("n_settings=%d", nrow(grid)))
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  step_id <- sprintf("grid_%02d", i)
  pars <- grid[i]
  if (i %% 5 == 0 || i == 1 || i == nrow(grid)) {
    log_msg(sprintf("Sweep %d/%d", i, nrow(grid)))
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
fwrite(res, file.path(OUT_TABLE, "qc_parameter_sweep_summary.tsv"), sep = "\t")

ok <- res[status == "ok"]
if (nrow(ok) > 0) {
  best <- ok[order(grey_pct, -non_grey_modules, -module_size_median)][1]
  fwrite(best, file.path(OUT_TABLE, "qc_parameter_sweep_recommended.tsv"), sep = "\t")
}

progress_update("global", "ok", sprintf("completed_n=%d", nrow(res)))
log_msg("Parameter sweep complete")
