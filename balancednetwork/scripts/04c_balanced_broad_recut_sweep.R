#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(callr)
  library(data.table)
  library(WGCNA)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)
allow_balanced_wgcna_threads()

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

setting_key <- function(power, deepSplit, mergeCutHeight, minModuleSize) {
  sprintf("%s|%s|%.12g|%s", power, deepSplit, mergeCutHeight, minModuleSize)
}

read_completed_keys <- function(path) {
  if (!file.exists(path)) return(character())
  existing <- fread(path)
  if (!nrow(existing)) return(character())
  existing[, setting_key(power, deepSplit, mergeCutHeight, minModuleSize)]
}

write_live_summary <- function(path, summary_path, total_settings) {
  if (!file.exists(path)) return(invisible(NULL))
  x <- fread(path)
  ok <- x[status == "ok"]
  candidates <- ok[
    non_grey_modules >= BAL_PARAMS$broad_keep_module_min &
      non_grey_modules <= BAL_PARAMS$broad_keep_module_max &
      grey_pct <= BAL_PARAMS$broad_keep_grey_pct_max &
      module_size_median >= BAL_PARAMS$broad_keep_module_median_min &
      largest_non_grey_module_pct <= BAL_PARAMS$broad_keep_largest_module_pct_max &
      singleton_non_grey_modules == 0
  ]
  best <- candidates[
    order(abs(non_grey_modules - BAL_PARAMS$broad_target_modules),
          grey_pct, -module_size_median, largest_non_grey_module_pct)
  ][1:min(10, .N)]

  summary <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    total_settings = total_settings,
    completed_settings = nrow(x),
    ok_settings = nrow(ok),
    failed_settings = x[status == "failed", .N],
    timeout_settings = x[status == "timeout", .N],
    candidate_settings = nrow(candidates)
  )

  fwrite(summary, summary_path, sep = "\t")
  if (nrow(best)) {
    fwrite(best, sub("\\.tsv$", "_best.tsv", summary_path), sep = "\t")
  }
  invisible(NULL)
}

run_recut_with_timeout <- function(mx, net, pars, timeout_sec, n_threads) {
  callr::r(
    func = function(mx, net, pars, n_threads) {
      suppressPackageStartupMessages(library(WGCNA))
      WGCNA::allowWGCNAThreads(nThreads = n_threads)
      fit <- WGCNA::recutConsensusTrees(
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
      )
      fit$colors
    },
    args = list(mx = mx, net = net, pars = as.list(pars), n_threads = n_threads),
    timeout = timeout_sec,
    show = FALSE
  )
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
summary_file <- file.path(BAL$qc_broad_dir, "broad_recut_live_summary.tsv")

timeout_sec <- suppressWarnings(as.integer(Sys.getenv("BALANCED_BROAD_RECUT_TIMEOUT_SEC", "600")))
if (!is.finite(timeout_sec) || timeout_sec < 1L) timeout_sec <- 600L
n_threads <- suppressWarnings(as.integer(Sys.getenv("BALANCED_WGCNA_THREADS", "8")))
if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 8L
force_run <- identical(Sys.getenv("BALANCED_BROAD_FORCE", unset = "0"), "1")

if (force_run) {
  if (file.exists(out_file)) file.remove(out_file)
  if (file.exists(cache_manifest_file)) file.remove(cache_manifest_file)
}

progress_update(
  "broad_recut",
  "start",
  sprintf("n_settings=%d timeout_sec=%d n_threads=%d force=%s", nrow(grid), timeout_sec, n_threads, force_run)
)

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
    completed_keys <- read_completed_keys(out_file)
    key <- setting_key(pars$power, pars$deepSplit, pars$mergeCutHeight, pars$minModuleSize)
    if (!force_run && key %in% completed_keys) {
      next
    }

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

    started_at <- Sys.time()
    colors <- tryCatch(
      run_recut_with_timeout(mx, net, pars, timeout_sec, n_threads),
      error = function(e) e
    )
    elapsed_sec <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

    if (inherits(colors, "error")) {
      err <- conditionMessage(colors)
      status <- if (grepl("timed out|timeout", err, ignore.case = TRUE)) "timeout" else "failed"
      row <- data.table(
        power = pars$power,
        deepSplit = pars$deepSplit,
        mergeCutHeight = pars$mergeCutHeight,
        minModuleSize = pars$minModuleSize,
        status = status,
        error = err,
        elapsed_sec = elapsed_sec
      )
      progress_update(recut_stage, status, sprintf("elapsed_sec=%.1f error=%s", elapsed_sec, err))
    } else {
      row <- cbind(
        data.table(
          power = pars$power,
          deepSplit = pars$deepSplit,
          mergeCutHeight = pars$mergeCutHeight,
          minModuleSize = pars$minModuleSize,
          status = "ok",
          error = "",
          elapsed_sec = elapsed_sec
        ),
        module_summary_from_colors(colors)
      )
      progress_update(
        recut_stage,
        "ok",
        sprintf(
          "non_grey=%d grey_pct=%.2f elapsed_sec=%.1f",
          row$non_grey_modules[[1]],
          row$grey_pct[[1]],
          elapsed_sec
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
    write_live_summary(out_file, summary_file, nrow(grid))
  }
}

progress_update("broad_recut", "ok", sprintf("saved=%d", nrow(grid)))
write_live_summary(out_file, summary_file, nrow(grid))
log_msg("balanced broad recut sweep complete")
