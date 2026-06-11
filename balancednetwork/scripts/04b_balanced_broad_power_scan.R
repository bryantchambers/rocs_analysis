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

ctx <- load_balanced_context()
expr_by_core <- ctx$expr_by_core
powers <- BAL_PARAMS$broad_power_scan

progress_update("power_scan", "start", sprintf("n_powers=%d", length(powers)))

scan_long <- rbindlist(lapply(names(expr_by_core), function(core_id) {
  progress_update("power_scan_core", "start", sprintf("core=%s", core_id))
  sft <- pickSoftThreshold(
    expr_by_core[[core_id]],
    powerVector = powers,
    networkType = "signed",
    verbose = 0
  )
  fit <- as.data.table(sft$fitIndices)
  setnames(fit, old = intersect(c("SFT.R.sq", "mean.k."), names(fit)),
           new = c(if ("SFT.R.sq" %in% names(fit)) "signedR2" else NULL,
                   if ("mean.k." %in% names(fit)) "meanK" else NULL),
           skip_absent = TRUE)
  fit[, core := core_id]
  progress_update("power_scan_core", "ok", sprintf("core=%s rows=%d", core_id, nrow(fit)))
  fit
}), fill = TRUE)

scan_summary <- scan_long[, .(
  mean_signedR2 = mean(signedR2, na.rm = TRUE),
  min_signedR2 = min(signedR2, na.rm = TRUE),
  mean_meanK = mean(meanK, na.rm = TRUE),
  min_meanK = min(meanK, na.rm = TRUE)
), by = Power][order(Power)]

scan_summary[, keep_for_recut := Power %in% BAL_PARAMS$broad_power_recut]
scan_summary[, keep_reason := fifelse(
  keep_for_recut,
  "in broad recut grid",
  "power scan only"
)]

dir.create(BAL$qc_broad_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(scan_long, file.path(BAL$qc_broad_dir, "broad_power_scan_long.tsv"), sep = "\t")
fwrite(scan_summary, file.path(BAL$qc_broad_dir, "broad_power_scan.tsv"), sep = "\t")

progress_update("power_scan", "ok", sprintf("saved_powers=%d", nrow(scan_summary)))
log_msg("balanced broad power scan complete")
