#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)
allow_balanced_wgcna_threads()

progress_tsv <- Sys.getenv("BALANCED_NEIGHBORHOOD_PROGRESS_TSV", unset = "")
if (!nzchar(progress_tsv)) {
  progress_tsv <- file.path(BAL$logs_dir, sprintf("neighborhood_scan_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
}

progress_update <- function(stage, status, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      lane = "balancednetwork_neighborhood_scan",
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

setting_key <- function(power, deepSplit, mergeCutHeight, minModuleSize) {
  sprintf("%s|%s|%.12g|%s", power, deepSplit, mergeCutHeight, minModuleSize)
}

read_completed_keys <- function(path) {
  if (!file.exists(path)) return(character())
  existing <- fread(path, fill = TRUE)
  if (!nrow(existing)) return(character())
  existing[, setting_key(power, deepSplit, mergeCutHeight, minModuleSize)]
}

make_cluster_grid <- function(cluster, purpose, power_vals, deep_vals, mms_vals) {
  CJ(
    cluster = cluster,
    purpose = purpose,
    power = as.integer(power_vals),
    deepSplit = as.integer(deep_vals),
    mergeCutHeight = BAL_PARAMS$neighborhood_mergeCutHeight,
    minModuleSize = as.integer(mms_vals),
    unique = TRUE
  )
}

cluster_grid <- rbindlist(list(
  make_cluster_grid("A", "primary 7-module low-power cluster", 2:4, 3:4, c(8L, 9L, 10L)),
  make_cluster_grid("B", "primary 7-module higher-power cluster", 7:9, 3:4, c(8L, 9L, 10L)),
  make_cluster_grid("C", "fallback 6-module family", 9:16, 2:4, c(8L, 9L, 10L)),
  make_cluster_grid("D", "upper-bound 8-module family", 4:6, 2:4, c(8L, 9L))
), fill = TRUE, use.names = TRUE)

grid <- cluster_grid[
  ,
  .(
    cluster_labels = paste(sort(unique(cluster)), collapse = ","),
    cluster_purposes = paste(unique(purpose), collapse = " | ")
  ),
  by = .(power, deepSplit, mergeCutHeight, minModuleSize)
]

limit_n <- suppressWarnings(as.integer(Sys.getenv("BALANCED_NEIGHBORHOOD_LIMIT", "")))
if (is.finite(limit_n) && limit_n > 0L && limit_n < nrow(grid)) {
  setorderv(grid, c("power", "deepSplit", "mergeCutHeight", "minModuleSize"))
  grid <- grid[1:limit_n]
}

ctx <- load_balanced_context()
expr_by_core <- ctx$expr_by_core

dir.create(BAL$qc_neighborhood_dir, recursive = TRUE, showWarnings = FALSE)
grid_file <- file.path(BAL$qc_neighborhood_dir, "neighborhood_grid.tsv")
out_file <- file.path(BAL$qc_neighborhood_dir, "neighborhood_scan.tsv")

force_run <- identical(Sys.getenv("BALANCED_NEIGHBORHOOD_FORCE", unset = "0"), "1")
if (force_run && file.exists(out_file)) file.remove(out_file)
fwrite(grid, grid_file, sep = "\t")

progress_update(
  "neighborhood_scan",
  "start",
  sprintf("n_settings=%d force=%s", nrow(grid), force_run)
)

for (i in seq_len(nrow(grid))) {
  pars <- grid[i]
  key <- setting_key(pars$power, pars$deepSplit, pars$mergeCutHeight, pars$minModuleSize)
  if (!force_run && key %in% read_completed_keys(out_file)) next

  stage_id <- sprintf("setting_%03d", i)
  progress_update(
    stage_id,
    "start",
    sprintf(
      "cluster=%s power=%d deepSplit=%d mergeCutHeight=%.2f minModuleSize=%d",
      pars$cluster_labels, pars$power, pars$deepSplit, pars$mergeCutHeight, pars$minModuleSize
    )
  )

  started_at <- Sys.time()
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
  elapsed_sec <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

  if (inherits(fit, "error")) {
    row <- data.table(
      cluster_labels = pars$cluster_labels,
      cluster_purposes = pars$cluster_purposes,
      power = pars$power,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      minModuleSize = pars$minModuleSize,
      status = "failed",
      error = conditionMessage(fit),
      elapsed_sec = elapsed_sec
    )
    progress_update(stage_id, "failed", sprintf("elapsed_sec=%.1f error=%s", elapsed_sec, conditionMessage(fit)))
  } else {
    row <- cbind(
      data.table(
        cluster_labels = pars$cluster_labels,
        cluster_purposes = pars$cluster_purposes,
        power = pars$power,
        deepSplit = pars$deepSplit,
        mergeCutHeight = pars$mergeCutHeight,
        minModuleSize = pars$minModuleSize,
        status = "ok",
        error = "",
        elapsed_sec = elapsed_sec
      ),
      module_summary_from_colors(fit$colors)
    )
    progress_update(
      stage_id,
      "ok",
      sprintf(
        "cluster=%s non_grey=%d grey_pct=%.2f elapsed_sec=%.1f",
        pars$cluster_labels,
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
}

progress_update("neighborhood_scan", "complete", sprintf("saved=%d", nrow(grid)))
log_msg("balanced neighborhood scan complete")
