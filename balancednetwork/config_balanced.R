#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

BAL <- list(
  base_dir = here("balancednetwork"),
  results_dir = here("balancednetwork", "results"),
  tables_dir = here("balancednetwork", "results", "tables"),
  wgcna_dir = here("balancednetwork", "results", "wgcna"),
  stability_dir = here("balancednetwork", "results", "stability"),
  reports_dir = here("balancednetwork", "results", "reports"),
  qc_dir = here("balancednetwork", "results", "qc"),
  qc_tables_dir = here("balancednetwork", "results", "qc", "tables"),
  qc_fig_dir = here("balancednetwork", "results", "qc", "figures"),
  qc_full_dir = here("balancednetwork", "results", "qc", "full_eval"),
  logs_dir = here("balancednetwork", "results", "logs")
)

for (d in BAL) dir.create(d, recursive = TRUE, showWarnings = FALSE)

BAL_PARAMS <- list(
  train_cores = PARAMS$stage1_cores,
  validation_core = PARAMS$validation_core,
  soft_power = 12L,
  deepSplit = 3L,
  mergeCutHeight = 0.25,
  minModuleSize = 20L,
  grid_power = c(12L, 16L, 20L),
  grid_deepSplit = c(1L, 2L, 3L),
  grid_mergeCutHeight = c(0.10, 0.15, 0.20, 0.25),
  grid_minModuleSize = c(20L, 30L),
  age_grid_points = PARAMS$wgcna_stability_age_grid_points
)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

norm01 <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(0.5, length(x)))
  z <- (x - rng[1]) / (rng[2] - rng[1])
  if (!higher_better) z <- 1 - z
  z
}

load_balanced_context <- function() {
  vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
  meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
  sel <- fread(file.path(BAL$tables_dir, "balanced_baseline_samples.tsv"))
  quotas <- fread(file.path(BAL$tables_dir, "balance_bin_quotas.tsv"))

  expr_by_core <- lapply(BAL_PARAMS$train_cores, function(core_id) {
    samps <- sel[core == core_id, label]
    samps <- intersect(samps, rownames(vst))
    vst[samps, , drop = FALSE]
  })
  names(expr_by_core) <- BAL_PARAMS$train_cores

  validation_samples <- intersect(meta[core == BAL_PARAMS$validation_core, label], rownames(vst))
  expr_validation <- vst[validation_samples, , drop = FALSE]

  pool_by_core_bin <- lapply(BAL_PARAMS$train_cores, function(core_id) {
    dt <- sel[core == core_id, .(label, age_bin)]
    split(dt$label, dt$age_bin)
  })
  names(pool_by_core_bin) <- BAL_PARAMS$train_cores

  list(
    vst = vst,
    meta = meta,
    sel = sel,
    quotas = quotas,
    expr_by_core = expr_by_core,
    expr_validation = expr_validation,
    pool_by_core_bin = pool_by_core_bin
  )
}

fit_consensus_from_expr <- function(expr_by_core, power, deepSplit, mergeCutHeight, minModuleSize) {
  mx <- lapply(BAL_PARAMS$train_cores, function(core_id) list(data = expr_by_core[[core_id]]))
  names(mx) <- BAL_PARAMS$train_cores
  blockwiseConsensusModules(
    multiExpr = mx,
    power = power,
    networkType = "signed",
    corType = "pearson",
    maxBlockSize = 5000,
    minModuleSize = minModuleSize,
    deepSplit = deepSplit,
    mergeCutHeight = mergeCutHeight,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )
}
