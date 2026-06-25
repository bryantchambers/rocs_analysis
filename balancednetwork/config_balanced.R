#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

BALANCED_BROAD_SCAN_RUN_ID <- Sys.getenv("BALANCED_BROAD_SCAN_RUN_ID", unset = "")
BALANCED_NEIGHBORHOOD_RUN_ID <- Sys.getenv("BALANCED_NEIGHBORHOOD_RUN_ID", unset = "")
BALANCED_STAGE1_VARIANT <- Sys.getenv("BALANCED_STAGE1_VARIANT", unset = "")
BALANCED_RUN_PROGRESS_TSV <- Sys.getenv("BALANCED_RUN_PROGRESS_TSV", unset = "")
BALANCED_NEIGHBORHOOD_EVAL_TOP_N <- suppressWarnings(as.integer(Sys.getenv("BALANCED_NEIGHBORHOOD_EVAL_TOP_N", unset = "8")))
if (!is.finite(BALANCED_NEIGHBORHOOD_EVAL_TOP_N) || BALANCED_NEIGHBORHOOD_EVAL_TOP_N < 1L) {
  BALANCED_NEIGHBORHOOD_EVAL_TOP_N <- 8L
}

if (nzchar(BALANCED_STAGE1_VARIANT)) {
  RESULTS$stage1 <- here("balancednetwork", "results", "stage1", BALANCED_STAGE1_VARIANT)
}

BAL <- list(
  base_dir = here("balancednetwork"),
  results_dir = here("balancednetwork", "results"),
  stage1_variant = if (nzchar(BALANCED_STAGE1_VARIANT)) BALANCED_STAGE1_VARIANT else "shared",
  stage1_dir = RESULTS$stage1,
  tables_dir = here("balancednetwork", "results", "tables"),
  wgcna_dir = here("balancednetwork", "results", "wgcna"),
  stability_dir = here("balancednetwork", "results", "stability"),
  reports_dir = here("balancednetwork", "results", "reports"),
  qc_dir = here("balancednetwork", "results", "qc"),
  qc_tables_dir = here("balancednetwork", "results", "qc", "tables"),
  qc_fig_dir = here("balancednetwork", "results", "qc", "figures"),
  qc_full_dir = here("balancednetwork", "results", "qc", "full_eval"),
  qc_neighborhood_dir = if (nzchar(BALANCED_NEIGHBORHOOD_RUN_ID)) {
    here("balancednetwork", "results", "qc", "neighborhood_scan", BALANCED_NEIGHBORHOOD_RUN_ID)
  } else {
    here("balancednetwork", "results", "qc", "neighborhood_scan")
  },
  qc_neighborhood_full_dir = if (nzchar(BALANCED_NEIGHBORHOOD_RUN_ID)) {
    here("balancednetwork", "results", "qc", "neighborhood_scan", BALANCED_NEIGHBORHOOD_RUN_ID, "full_eval")
  } else {
    here("balancednetwork", "results", "qc", "neighborhood_scan", "full_eval")
  },
  qc_broad_dir = if (nzchar(BALANCED_BROAD_SCAN_RUN_ID)) {
    here("balancednetwork", "results", "qc", "broad_sweep", BALANCED_BROAD_SCAN_RUN_ID)
  } else {
    here("balancednetwork", "results", "qc", "broad_sweep")
  },
  qc_broad_cache_dir = if (nzchar(BALANCED_BROAD_SCAN_RUN_ID)) {
    here("balancednetwork", "results", "qc", "broad_sweep", BALANCED_BROAD_SCAN_RUN_ID, "cache")
  } else {
    here("balancednetwork", "results", "qc", "broad_sweep", "cache")
  },
  logs_dir = here("balancednetwork", "results", "logs")
)

for (d in BAL[c(
  "base_dir",
  "results_dir",
  "stage1_dir",
  "tables_dir",
  "wgcna_dir",
  "stability_dir",
  "reports_dir",
  "qc_dir",
  "qc_tables_dir",
  "qc_fig_dir",
  "qc_full_dir",
  "qc_neighborhood_dir",
  "qc_neighborhood_full_dir",
  "qc_broad_dir",
  "qc_broad_cache_dir",
  "logs_dir"
)]) dir.create(d, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS$stage1, recursive = TRUE, showWarnings = FALSE)

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
  broad_power_scan = c(1:10, 12L, 14L, 16L),
  broad_power_recut = c(2:10, 12L, 14L, 16L),
  broad_deepSplit = 0:4,
  broad_mergeCutHeight = c(0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20),
  broad_minModuleSize = c(8L, 10L, 12L, 15L, 20L),
  broad_target_modules = 6L,
  broad_keep_module_min = 4L,
  broad_keep_module_max = 9L,
  broad_keep_grey_pct_max = 70,
  broad_keep_module_median_min = 10,
  broad_keep_largest_module_pct_max = 45,
  broad_top_n_confirm = 15L,
  neighborhood_mergeCutHeight = c(0.01, 0.02, 0.05, 0.08),
  neighborhood_main_grey_pct_max = 60,
  neighborhood_secondary_grey_pct_max = 65,
  neighborhood_largest_module_pct_max = 20,
  neighborhood_eval_top_n = BALANCED_NEIGHBORHOOD_EVAL_TOP_N,
  age_grid_points = PARAMS$wgcna_stability_age_grid_points
)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

allow_balanced_wgcna_threads <- function(default_threads = 8L) {
  n_threads <- suppressWarnings(as.integer(Sys.getenv("BALANCED_WGCNA_THREADS", default_threads)))
  if (!is.finite(n_threads) || n_threads < 1L) n_threads <- default_threads
  WGCNA::allowWGCNAThreads(nThreads = n_threads)
  invisible(n_threads)
}

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

module_summary_from_colors <- function(colors) {
  tab <- sort(table(as.character(colors)), decreasing = TRUE)
  non_g_names <- setdiff(names(tab), c("grey", "gold"))
  non_g_sizes <- as.numeric(tab[non_g_names])
  grey_n <- if ("grey" %in% names(tab)) as.numeric(tab[["grey"]]) else 0
  total_n <- length(colors)
  data.table(
    non_grey_modules = length(non_g_names),
    grey_pct = grey_n / total_n * 100,
    module_size_median = if (length(non_g_sizes)) median(non_g_sizes) else NA_real_,
    min_non_grey_module_size = if (length(non_g_sizes)) min(non_g_sizes) else NA_real_,
    max_non_grey_module_size = if (length(non_g_sizes)) max(non_g_sizes) else NA_real_,
    largest_non_grey_module_pct = if (length(non_g_sizes)) max(non_g_sizes) / total_n * 100 else NA_real_,
    singleton_non_grey_modules = sum(non_g_sizes <= 1),
    total_taxa = total_n
  )
}
