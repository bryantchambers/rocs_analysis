#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)

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

rank_bucket <- function(non_grey_modules) {
  fifelse(
    non_grey_modules == 7L, 0L,
    fifelse(non_grey_modules == 6L, 1L,
            fifelse(non_grey_modules == 5L, 2L,
                    fifelse(non_grey_modules == 8L, 3L, 4L)))
  )
}

preferred_merge_rank <- function(x) {
  vals <- c(0.02, 0.05, 0.08, 0.01)
  idx <- match(round(as.numeric(x), 2), vals)
  idx[is.na(idx)] <- length(vals) + 1L
  idx
}

sweep_file <- file.path(BAL$qc_neighborhood_dir, "neighborhood_scan.tsv")
if (!file.exists(sweep_file)) stop("Missing neighborhood scan: ", sweep_file)
progress_update("neighborhood_candidates", "start")

sweep <- fread(sweep_file, fill = TRUE)
if (nrow(sweep) == 0L) stop("Neighborhood scan is empty: ", sweep_file)

sweep[, raw_row_index := .I]
setorderv(sweep, c("power", "deepSplit", "mergeCutHeight", "minModuleSize", "raw_row_index"))
sweep_latest <- sweep[, .SD[.N], by = .(power, deepSplit, mergeCutHeight, minModuleSize)]
ok <- copy(sweep_latest[status == "ok"])
if (nrow(ok) == 0L) stop("No successful neighborhood scan rows.")

ok[, cluster_primary := tstrsplit(cluster_labels, ",", fixed = TRUE, keep = 1L)]
ok[, main_candidate := (
  non_grey_modules >= 5L &
    non_grey_modules <= 7L &
    minModuleSize %in% c(8L, 9L, 10L) &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct <= BAL_PARAMS$neighborhood_largest_module_pct_max &
    grey_pct <= BAL_PARAMS$neighborhood_main_grey_pct_max
)]
ok[, secondary_candidate := (
  non_grey_modules >= 5L &
    non_grey_modules <= 7L &
    minModuleSize %in% c(8L, 9L, 10L) &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct <= BAL_PARAMS$neighborhood_largest_module_pct_max &
    grey_pct <= BAL_PARAMS$neighborhood_secondary_grey_pct_max
)]
ok[, boundary_candidate := (
  non_grey_modules == 8L &
    minModuleSize %in% c(8L, 9L) &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct <= BAL_PARAMS$neighborhood_largest_module_pct_max
)]

ok[, merge_phenotype_key := paste(
  non_grey_modules,
  sprintf("%.2f", round(grey_pct, 2)),
  sprintf("%.1f", round(module_size_median, 1)),
  sprintf("%.2f", round(largest_non_grey_module_pct, 2)),
  sep = "|"
)]

regions <- ok[
  ,
  .(
    cluster_labels = first(cluster_labels),
    cluster_primary = first(cluster_primary),
    cluster_purposes = first(cluster_purposes),
    merge_support_n = uniqueN(mergeCutHeight),
    mergeCutHeight_values = paste(sprintf("%.2f", sort(unique(mergeCutHeight))), collapse = ","),
    mergeCutHeight_min = min(mergeCutHeight),
    mergeCutHeight_max = max(mergeCutHeight),
    representative_mergeCutHeight = mergeCutHeight[order(preferred_merge_rank(mergeCutHeight), mergeCutHeight)][1],
    non_grey_modules = first(non_grey_modules),
    grey_pct = first(grey_pct),
    module_size_median = first(module_size_median),
    min_non_grey_module_size = first(min_non_grey_module_size),
    max_non_grey_module_size = first(max_non_grey_module_size),
    largest_non_grey_module_pct = first(largest_non_grey_module_pct),
    singleton_non_grey_modules = first(singleton_non_grey_modules),
    total_taxa = first(total_taxa),
    elapsed_sec_median = median(elapsed_sec, na.rm = TRUE),
    main_candidate = all(main_candidate),
    secondary_candidate = all(secondary_candidate),
    boundary_candidate = all(boundary_candidate)
  ),
  by = .(power, deepSplit, minModuleSize, merge_phenotype_key)
]

regions[, module_rank_bucket := rank_bucket(non_grey_modules)]
regions[, cluster_priority := match(cluster_primary, c("A", "B", "C", "D"))]
regions[is.na(cluster_priority), cluster_priority := 99L]

setorder(
  regions,
  -main_candidate,
  -secondary_candidate,
  -boundary_candidate,
  module_rank_bucket,
  grey_pct,
  -merge_support_n,
  -module_size_median,
  largest_non_grey_module_pct,
  cluster_priority,
  power,
  deepSplit,
  minModuleSize,
  representative_mergeCutHeight
)
regions[, candidate_rank := .I]
regions[, setting_id := sprintf("nh_%03d", .I)]

quotas <- data.table(cluster_primary = c("A", "B", "C", "D"), quota = c(2L, 2L, 2L, 1L))
eligible_regions <- regions[main_candidate == TRUE | secondary_candidate == TRUE | boundary_candidate == TRUE]
selection_pool <- if (nrow(eligible_regions) > 0L) eligible_regions else regions

selected_by_cluster <- rbindlist(lapply(seq_len(nrow(quotas)), function(i) {
  q <- quotas[i]
  pool <- selection_pool[cluster_primary == q$cluster_primary]
  if (nrow(pool) == 0L) return(selection_pool[0])
  pool[1:min(q$quota, .N)]
}), fill = TRUE)

remaining <- selection_pool[!setting_id %in% selected_by_cluster$setting_id]
fill_n <- max(0L, BAL_PARAMS$neighborhood_eval_top_n - nrow(selected_by_cluster))
selected_fill <- if (fill_n > 0L && nrow(remaining) > 0L) remaining[1:min(fill_n, .N)] else remaining[0]
settings_to_eval <- rbindlist(
  Filter(function(x) is.data.table(x) && nrow(x) > 0L, list(selected_by_cluster, selected_fill)),
  fill = TRUE
)
settings_to_eval <- unique(settings_to_eval, by = "setting_id")
settings_to_eval[, mergeCutHeight := representative_mergeCutHeight]
setorder(settings_to_eval, candidate_rank)
if (is.finite(BAL_PARAMS$neighborhood_eval_top_n) && BAL_PARAMS$neighborhood_eval_top_n > 0L && nrow(settings_to_eval) > BAL_PARAMS$neighborhood_eval_top_n) {
  settings_to_eval <- settings_to_eval[seq_len(BAL_PARAMS$neighborhood_eval_top_n)]
}
if ("quota" %in% names(settings_to_eval)) settings_to_eval[, quota := NULL]

candidate_rows <- regions[
  ok,
  on = .(
    power,
    deepSplit,
    minModuleSize,
    merge_phenotype_key,
    representative_mergeCutHeight = mergeCutHeight
  )
]

fwrite(regions, file.path(BAL$qc_neighborhood_dir, "neighborhood_candidates.tsv"), sep = "\t")
fwrite(settings_to_eval, file.path(BAL$qc_neighborhood_dir, "neighborhood_settings_to_full_eval.tsv"), sep = "\t")
fwrite(candidate_rows, file.path(BAL$qc_neighborhood_dir, "neighborhood_representative_rows.tsv"), sep = "\t")

report_file <- file.path(BAL$qc_neighborhood_dir, "BALANCED_NEIGHBORHOOD_SCAN_REPORT.md")
sink(report_file)
cat("# Balanced Neighborhood Scan Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
if (nzchar(BALANCED_NEIGHBORHOOD_RUN_ID)) {
  cat("- Neighborhood run ID: `", BALANCED_NEIGHBORHOOD_RUN_ID, "`\n", sep = "")
}
if (nzchar(BALANCED_BROAD_SCAN_RUN_ID)) {
  cat("- Broad-scan seed source: `", BALANCED_BROAD_SCAN_RUN_ID, "`\n", sep = "")
}
cat("\n")
cat("## Scan Summary\n\n")
cat(sprintf("- Unique scanned settings: %d\n", nrow(sweep_latest)))
cat(sprintf("- Successful settings: %d\n", nrow(ok)))
cat(sprintf("- Failed settings: %d\n", sweep_latest[status == "failed", .N]))
cat(sprintf("- Main candidates: %d\n", regions[main_candidate == TRUE, .N]))
cat(sprintf("- Secondary candidates: %d\n", regions[secondary_candidate == TRUE, .N]))
cat(sprintf("- Boundary candidates: %d\n\n", regions[boundary_candidate == TRUE, .N]))

cat("## Top Neighborhood Candidates\n\n")
cat("|setting_id|cluster|power|deepSplit|mergeCutHeight|minModuleSize|non_grey_modules|grey_pct|module_size_median|largest_non_grey_module_pct|merge_support_n|\n")
cat("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
show_dt <- regions[1:min(12, .N)]
for (i in seq_len(nrow(show_dt))) {
  r <- show_dt[i]
  cat(sprintf(
    "|%s|%s|%d|%d|%.2f|%d|%d|%.2f|%.1f|%.2f|%d|\n",
    r$setting_id,
    r$cluster_labels,
    r$power,
    r$deepSplit,
    r$representative_mergeCutHeight,
    r$minModuleSize,
    r$non_grey_modules,
    r$grey_pct,
    r$module_size_median,
    r$largest_non_grey_module_pct,
    r$merge_support_n
  ))
}
cat("\n")

cat("## Settings Selected For Full Evaluation\n\n")
cat("|setting_id|cluster|power|deepSplit|mergeCutHeight|minModuleSize|non_grey_modules|grey_pct|\n")
cat("|---|---|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(settings_to_eval))) {
  r <- settings_to_eval[i]
  cat(sprintf(
    "|%s|%s|%d|%d|%.2f|%d|%d|%.2f|\n",
    r$setting_id, r$cluster_labels, r$power, r$deepSplit,
    r$representative_mergeCutHeight, r$minModuleSize,
    r$non_grey_modules, r$grey_pct
  ))
}
sink()

progress_update("neighborhood_candidates", "complete", sprintf("selected_n=%d", nrow(settings_to_eval)))
log_msg("balanced neighborhood candidate ranking complete")
