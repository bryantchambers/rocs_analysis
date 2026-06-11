#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()

progress_tsv <- Sys.getenv("BALANCED_BROAD_SWEEP_PROGRESS_TSV", unset = "")
if (!nzchar(progress_tsv)) {
  progress_tsv <- file.path(BAL$logs_dir, sprintf("broad_sweep_progress_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
}

progress_update <- function(stage, status, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      lane = "balancednetwork_broad_sweep",
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

mx_from_expr <- function(expr_by_core) {
  mx <- lapply(BAL_PARAMS$train_cores, function(core_id) list(data = expr_by_core[[core_id]]))
  names(mx) <- BAL_PARAMS$train_cores
  mx
}

ctx <- load_balanced_context()
expr_by_core <- ctx$expr_by_core
mx <- mx_from_expr(expr_by_core)

grid <- CJ(
  power = BAL_PARAMS$broad_power_recut,
  deepSplit = BAL_PARAMS$broad_deepSplit,
  mergeCutHeight = BAL_PARAMS$broad_mergeCutHeight,
  minModuleSize = BAL_PARAMS$broad_minModuleSize
)

dir.create(BAL$qc_broad_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(BAL$qc_broad_cache_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(BAL$qc_broad_dir, "broad_recut_sweep.tsv")
cache_manifest_file <- file.path(BAL$qc_broad_dir, "broad_recut_cache_manifest.tsv")
if (file.exists(out_file)) file.remove(out_file)
if (file.exists(cache_manifest_file)) file.remove(cache_manifest_file)

progress_update("broad_recut", "start", sprintf("n_settings=%d", nrow(grid)))

for (pwr in BAL_PARAMS$broad_power_recut) {
  power_stage <- sprintf("power_%02d_cache", pwr)
  power_dir <- file.path(BAL$qc_broad_cache_dir, sprintf("power_%02d", pwr))
  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)

  progress_update(power_stage, "start", sprintf("power=%d", pwr))
  net <- blockwiseConsensusModules(
    multiExpr = mx,
    power = pwr,
    networkType = "signed",
    corType = "pearson",
    maxBlockSize = 5000,
    minModuleSize = min(BAL_PARAMS$broad_minModuleSize),
    deepSplit = 2,
    mergeCutHeight = 0.10,
    numericLabels = FALSE,
    saveIndividualTOMs = TRUE,
    saveConsensusTOMs = TRUE,
    cacheDir = power_dir,
    individualTOMFileNames = file.path(power_dir, "individualTOM-Set%s-Block%b.RData"),
    consensusTOMFilePattern = file.path(power_dir, "consensusTOM-block.%b.RData"),
    verbose = 0
  )

  cache_row <- data.table(
    power = pwr,
    consensus_tom_files = paste(net$consensusTOMInfo$TOMFiles, collapse = ";"),
    n_blocks = length(net$dendrograms),
    n_taxa = sum(net$goodGenes)
  )
  fwrite(
    cache_row,
    cache_manifest_file,
    sep = "\t",
    append = file.exists(cache_manifest_file),
    col.names = !file.exists(cache_manifest_file)
  )
  progress_update(power_stage, "ok", sprintf("power=%d blocks=%d", pwr, length(net$dendrograms)))

  power_grid <- grid[power == pwr]
  for (i in seq_len(nrow(power_grid))) {
    pars <- power_grid[i]
    recut_stage <- sprintf(
      "power_%02d_recut_%03d",
      pwr,
      i
    )
    progress_update(
      recut_stage,
      "start",
      sprintf(
        "power=%d deepSplit=%d mergeCutHeight=%.2f minModuleSize=%d",
        pars$power, pars$deepSplit, pars$mergeCutHeight, pars$minModuleSize
      )
    )
    fit <- tryCatch(
      recutConsensusTrees(
        multiExpr = mx,
        goodSamples = net$goodSamples,
        goodGenes = net$goodGenes,
        blocks = net$blocks,
        TOMFiles = net$consensusTOMInfo$TOMFiles,
        dendrograms = net$dendrograms,
        corType = "pearson",
        networkType = "signed",
        deepSplit = pars$deepSplit,
        minModuleSize = pars$minModuleSize,
        mergeCutHeight = pars$mergeCutHeight,
        numericLabels = FALSE,
        verbose = 0
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      row <- data.table(
        power = pars$power,
        deepSplit = pars$deepSplit,
        mergeCutHeight = pars$mergeCutHeight,
        minModuleSize = pars$minModuleSize,
        status = "failed",
        error = conditionMessage(fit)
      )
      progress_update(recut_stage, "failed", conditionMessage(fit))
    } else {
      row <- cbind(
        data.table(
          power = pars$power,
          deepSplit = pars$deepSplit,
          mergeCutHeight = pars$mergeCutHeight,
          minModuleSize = pars$minModuleSize,
          status = "ok",
          error = ""
        ),
        module_summary_from_colors(fit$colors)
      )
      progress_update(
        recut_stage,
        "ok",
        sprintf(
          "non_grey=%d grey_pct=%.2f",
          row$non_grey_modules[[1]],
          row$grey_pct[[1]]
        )
      )
    }

    fwrite(
      row,
      out_file,
      sep = "\t",
      append = file.exists(out_file),
      col.names = !file.exists(out_file)
    )
  }
}

progress_update("broad_recut", "ok", sprintf("saved=%d", nrow(grid)))
log_msg("balanced broad recut sweep complete")
