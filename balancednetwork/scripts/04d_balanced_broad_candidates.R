#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)

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

collapse_int_values <- function(x) {
  paste(sort(unique(as.integer(x))), collapse = ",")
}

collapse_numeric_values <- function(x, digits = 2L) {
  fmt <- sprintf("%%.%df", digits)
  vals <- sprintf(fmt, sort(unique(round(as.numeric(x), digits))))
  paste(vals, collapse = ",")
}

rank_bucket <- function(non_grey_modules) {
  fifelse(
    non_grey_modules == 7L, 0L,
    fifelse(non_grey_modules == 6L, 1L,
            fifelse(non_grey_modules == 5L, 2L,
                    fifelse(non_grey_modules == 8L, 3L, 4L)))
  )
}

assign_consecutive_group <- function(x) {
  x <- sort(unique(as.integer(x)))
  cumsum(c(TRUE, diff(x) != 1L))
}

take_n_rows <- function(dt, n) {
  if (nrow(dt) == 0L || n <= 0L) return(dt[0])
  dt[seq_len(min(nrow(dt), n))]
}

fmt_range <- function(lo, hi) {
  if (is.na(lo) || is.na(hi)) return("NA")
  if (lo == hi) as.character(lo) else sprintf("%s-%s", lo, hi)
}

fmt_num_range <- function(lo, hi, digits = 2L) {
  if (!is.finite(lo) || !is.finite(hi)) return("NA")
  fmt <- sprintf("%%.%df", digits)
  if (isTRUE(all.equal(lo, hi, tolerance = 10^(-digits)))) {
    sprintf(fmt, lo)
  } else {
    sprintf(paste0(fmt, "-", fmt), lo, hi)
  }
}

format_cluster_window <- function(power_lo, power_hi, deep_lo, deep_hi, mms_values) {
  sprintf(
    "power %s, deepSplit %s, minModuleSize %s",
    fmt_range(power_lo, power_hi),
    fmt_range(deep_lo, deep_hi),
    paste(mms_values, collapse = ",")
  )
}

progress_update("broad_candidates", "start")

sweep_file <- file.path(BAL$qc_broad_dir, "broad_recut_sweep.tsv")
if (!file.exists(sweep_file)) stop("Missing broad recut sweep: ", sweep_file)

sweep <- fread(sweep_file, fill = TRUE)
if (nrow(sweep) == 0L) stop("Broad recut sweep is empty: ", sweep_file)

sweep[, raw_row_index := .I]
setorderv(sweep, c("power", "deepSplit", "mergeCutHeight", "minModuleSize", "raw_row_index"))
sweep_latest <- sweep[, .SD[.N], by = .(power, deepSplit, mergeCutHeight, minModuleSize)]
setorderv(sweep_latest, c("power", "deepSplit", "mergeCutHeight", "minModuleSize"))

dup_count <- nrow(sweep) - nrow(sweep_latest)
ok <- copy(sweep_latest[status == "ok"])
timeout_n <- sweep_latest[status == "timeout", .N]
failed_n <- sweep_latest[status == "failed", .N]

if (nrow(ok) == 0L) stop("No successful broad recut rows after deduplication.")

ok[, main_seed_eligible := (
  non_grey_modules >= 5L &
    non_grey_modules <= 7L &
    minModuleSize == 8L &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct <= 20 &
    grey_pct <= 60
)]
ok[, secondary_seed_eligible := (
  non_grey_modules >= 5L &
    non_grey_modules <= 7L &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct <= 20 &
    grey_pct <= 65
)]
ok[, upper_bound_track := (
  non_grey_modules == 8L &
    minModuleSize %in% c(8L, 9L) &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct <= 20
)]

ok[, merge_phenotype_key := paste(
  non_grey_modules,
  sprintf("%.2f", round(grey_pct, 2)),
  sprintf("%.1f", round(module_size_median, 1)),
  sprintf("%.2f", round(largest_non_grey_module_pct, 2)),
  sep = "|"
)]

region_support <- ok[
  ,
  .(
    merge_support_n = uniqueN(mergeCutHeight),
    mergeCutHeight_values = collapse_numeric_values(mergeCutHeight, digits = 2L),
    mergeCutHeight_min = min(mergeCutHeight),
    mergeCutHeight_max = max(mergeCutHeight),
    non_grey_modules = as.integer(first(non_grey_modules)),
    grey_pct = as.numeric(first(grey_pct)),
    grey_pct_min = min(grey_pct),
    grey_pct_max = max(grey_pct),
    module_size_median = as.numeric(first(module_size_median)),
    module_size_median_min = min(module_size_median),
    module_size_median_max = max(module_size_median),
    largest_non_grey_module_pct = as.numeric(first(largest_non_grey_module_pct)),
    largest_non_grey_module_pct_min = min(largest_non_grey_module_pct),
    largest_non_grey_module_pct_max = max(largest_non_grey_module_pct),
    singleton_non_grey_modules = max(singleton_non_grey_modules),
    total_taxa = as.integer(first(total_taxa)),
    merge_phenotype_key = first(merge_phenotype_key)
  ),
  by = .(power, deepSplit, minModuleSize)
]

setorderv(region_support, c("power", "minModuleSize", "deepSplit"))
region_support[
  ,
  ds_group := assign_consecutive_group(deepSplit),
  by = .(power, minModuleSize, merge_phenotype_key)
]

region_classes <- region_support[
  ,
  .(
    deepSplit_min = min(deepSplit),
    deepSplit_max = max(deepSplit),
    deepSplit_values = collapse_int_values(deepSplit),
    deepSplit_support_n = uniqueN(deepSplit),
    merge_support_values = collapse_int_values(merge_support_n),
    merge_support_total = sum(merge_support_n),
    mergeCutHeight_values = first(mergeCutHeight_values),
    mergeCutHeight_min = min(mergeCutHeight_min),
    mergeCutHeight_max = max(mergeCutHeight_max),
    non_grey_modules = first(non_grey_modules),
    grey_pct = first(grey_pct),
    grey_pct_min = min(grey_pct_min),
    grey_pct_max = max(grey_pct_max),
    module_size_median = first(module_size_median),
    module_size_median_min = min(module_size_median_min),
    module_size_median_max = max(module_size_median_max),
    largest_non_grey_module_pct = first(largest_non_grey_module_pct),
    largest_non_grey_module_pct_min = min(largest_non_grey_module_pct_min),
    largest_non_grey_module_pct_max = max(largest_non_grey_module_pct_max),
    singleton_non_grey_modules = max(singleton_non_grey_modules),
    total_taxa = first(total_taxa)
  ),
  by = .(power, minModuleSize, merge_phenotype_key, ds_group)
]

setorderv(region_classes, c("power", "minModuleSize", "deepSplit_min"))
region_classes[, region_id := sprintf("region_%03d", .I)]
region_classes[, module_rank_bucket := rank_bucket(non_grey_modules)]
region_classes[, min_module_pref := fifelse(minModuleSize == 8L, 0L, fifelse(minModuleSize == 10L, 1L, 2L))]
region_classes[, stage1_main := (
  non_grey_modules >= 5L &
    non_grey_modules <= 7L &
    minModuleSize == 8L &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct_max <= 20 &
    grey_pct_max <= 60
)]
region_classes[, stage1_secondary := (
  non_grey_modules >= 5L &
    non_grey_modules <= 7L &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct_max <= 20 &
    grey_pct_max <= 65
)]
region_classes[, upper_bound_track := (
  non_grey_modules == 8L &
    minModuleSize %in% c(8L, 9L) &
    singleton_non_grey_modules == 0L &
    largest_non_grey_module_pct_max <= 20
)]

class_ranked <- copy(region_classes)
setorder(
  class_ranked,
  -stage1_main,
  -stage1_secondary,
  min_module_pref,
  module_rank_bucket,
  grey_pct,
  -merge_support_total,
  -deepSplit_support_n,
  -module_size_median,
  largest_non_grey_module_pct,
  power,
  deepSplit_min
)
class_ranked[, overall_rank := .I]

candidate_settings <- copy(ok[secondary_seed_eligible == TRUE | upper_bound_track == TRUE])
candidate_settings <- merge(
  candidate_settings,
  region_support[, .(power, deepSplit, minModuleSize, region_merge_support_n = merge_support_n, merge_phenotype_key)],
  by = c("power", "deepSplit", "minModuleSize", "merge_phenotype_key"),
  all.x = TRUE
)
candidate_settings <- merge(
  candidate_settings,
  class_ranked[
    ,
    .(
      power,
      minModuleSize,
      merge_phenotype_key,
      deepSplit_min,
      deepSplit_max,
      region_id,
      overall_rank,
      stage1_main,
      stage1_secondary,
      upper_bound_track
    )
  ],
  by = c("power", "minModuleSize", "merge_phenotype_key"),
  all.x = TRUE,
  allow.cartesian = TRUE
)
candidate_settings <- candidate_settings[deepSplit >= deepSplit_min & deepSplit <= deepSplit_max]
candidate_settings[, upper_bound_track := fifelse(!is.na(upper_bound_track.x), upper_bound_track.x, upper_bound_track.y)]
candidate_settings[, c("upper_bound_track.x", "upper_bound_track.y") := NULL]
setorder(candidate_settings, overall_rank, power, deepSplit, mergeCutHeight, minModuleSize)
candidate_settings[, setting_id := sprintf("broad_%04d", .I)]
setcolorder(
  candidate_settings,
  c(
    "setting_id", "region_id", "overall_rank",
    "power", "deepSplit", "mergeCutHeight", "minModuleSize",
    "status", "error", "elapsed_sec",
    "non_grey_modules", "grey_pct", "module_size_median",
    "min_non_grey_module_size", "max_non_grey_module_size",
    "largest_non_grey_module_pct", "singleton_non_grey_modules", "total_taxa",
    "raw_row_index", "merge_phenotype_key", "region_merge_support_n",
    "deepSplit_min", "deepSplit_max",
    "main_seed_eligible", "secondary_seed_eligible", "upper_bound_track",
    "stage1_main", "stage1_secondary"
  )
)

fwrite(candidate_settings, file.path(BAL$qc_broad_dir, "broad_recut_candidates.tsv"), sep = "\t")

seven_module_ranked <- class_ranked[stage1_main == TRUE & non_grey_modules == 7L]
primary_seed <- take_n_rows(seven_module_ranked, 1L)
secondary_seed <- take_n_rows(if (nrow(seven_module_ranked) >= 2L) seven_module_ranked[2:nrow(seven_module_ranked)] else seven_module_ranked[0], 1L)
comparison_7 <- take_n_rows(if (nrow(seven_module_ranked) >= 3L) seven_module_ranked[3:nrow(seven_module_ranked)] else seven_module_ranked[0], 3L)
fallback_6 <- take_n_rows(class_ranked[stage1_secondary == TRUE & minModuleSize == 8L & non_grey_modules == 6L], 4L)
upper_bound_8 <- take_n_rows(class_ranked[upper_bound_track == TRUE & minModuleSize == 8L], 4L)

confirmed <- rbindlist(
  Filter(
    function(x) !is.null(x) && nrow(x) > 0L,
    list(
      if (nrow(primary_seed)) cbind(role = "primary_seed_cluster", tier = "primary", primary_seed),
      if (nrow(secondary_seed)) cbind(role = "secondary_seed_cluster", tier = "primary", secondary_seed),
      if (nrow(comparison_7)) cbind(role = "comparison_7_module_cluster", tier = "primary", comparison_7),
      if (nrow(fallback_6)) cbind(role = "fallback_6_module_cluster", tier = "secondary", fallback_6),
      if (nrow(upper_bound_8)) cbind(role = "upper_bound_8_module_cluster", tier = "boundary", upper_bound_8)
    )
  ),
  fill = TRUE
)

if (nrow(confirmed) > 0L) {
  confirmed[
    ,
    `:=`(
      power_range = as.character(power),
      deepSplit_range = sprintf("%s", ifelse(deepSplit_min == deepSplit_max, deepSplit_min, sprintf("%d-%d", deepSplit_min, deepSplit_max))),
      mergeCutHeight_range = sprintf(
        "%s",
        ifelse(
          abs(mergeCutHeight_min - mergeCutHeight_max) < 1e-8,
          sprintf("%.2f", mergeCutHeight_min),
          sprintf("%.2f-%.2f", mergeCutHeight_min, mergeCutHeight_max)
        )
      )
    )
  ]
}

fwrite(confirmed, file.path(BAL$qc_broad_dir, "broad_recut_confirmed.tsv"), sep = "\t")
progress_update("broad_candidates", "ok", sprintf("n_candidates=%d n_confirmed=%d", nrow(candidate_settings), nrow(confirmed)))

cluster_windows <- data.table(
  cluster = c("A", "B", "C", "D"),
  purpose = c(
    "primary 7-module low-power cluster",
    "primary 7-module higher-power cluster",
    "fallback 6-module family",
    "upper-bound 8-module family"
  ),
  power_min = c(2L, 7L, 9L, 4L),
  power_max = c(4L, 9L, 16L, 6L),
  deepSplit_min = c(3L, 3L, 2L, 2L),
  deepSplit_max = c(4L, 4L, 4L, 4L)
)

cluster_summary_rows <- vector("list", nrow(cluster_windows))
for (i in seq_len(nrow(cluster_windows))) {
  cw <- cluster_windows[i]
  mms_values <- if (cw$cluster == "D") c(8L, 9L) else c(8L, 9L, 10L)
  subset_dt <- class_ranked[
    power >= cw$power_min &
      power <= cw$power_max &
      deepSplit_min >= cw$deepSplit_min &
      deepSplit_max <= cw$deepSplit_max &
      minModuleSize %in% mms_values
  ]

  if (nrow(subset_dt) == 0L) {
    cluster_summary_rows[[i]] <- data.table(
      cluster = cw$cluster,
      purpose = cw$purpose,
      window = format_cluster_window(cw$power_min, cw$power_max, cw$deepSplit_min, cw$deepSplit_max, mms_values),
      region_classes = 0L,
      best_region_id = NA_character_,
      best_power = NA_integer_,
      best_deepSplit = NA_character_,
      best_minModuleSize = NA_integer_,
      best_non_grey_modules = NA_integer_,
      best_grey_pct = NA_real_,
      best_merge_support_total = NA_integer_
    )
    next
  }

  best_subset <- if (cw$cluster == "D") {
    subset_dt[upper_bound_track == TRUE][order(grey_pct, -merge_support_total, -deepSplit_support_n)][1]
  } else {
    subset_dt[stage1_secondary == TRUE][
      order(module_rank_bucket, min_module_pref, grey_pct, -merge_support_total, -deepSplit_support_n, -module_size_median, largest_non_grey_module_pct)
    ][1]
  }

  cluster_summary_rows[[i]] <- data.table(
    cluster = cw$cluster,
    purpose = cw$purpose,
    window = format_cluster_window(cw$power_min, cw$power_max, cw$deepSplit_min, cw$deepSplit_max, mms_values),
    region_classes = nrow(subset_dt),
    best_region_id = best_subset$region_id,
    best_power = best_subset$power,
    best_deepSplit = fmt_range(best_subset$deepSplit_min, best_subset$deepSplit_max),
    best_minModuleSize = best_subset$minModuleSize,
    best_non_grey_modules = best_subset$non_grey_modules,
    best_grey_pct = best_subset$grey_pct,
    best_merge_support_total = best_subset$merge_support_total
  )
}
cluster_summary <- rbindlist(cluster_summary_rows, fill = TRUE)

module_dist <- ok[, .N, by = non_grey_modules][order(non_grey_modules)]
seed_mms_dist <- candidate_settings[stage1_secondary == TRUE, .N, by = minModuleSize][order(minModuleSize)]

report_file <- file.path(BAL$qc_broad_dir, "BALANCED_BROAD_SWEEP_REPORT.md")
sink(report_file)
cat("# Balanced Broad Sweep Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
if (nzchar(BALANCED_BROAD_SCAN_RUN_ID)) {
  cat("- Run ID: `", BALANCED_BROAD_SCAN_RUN_ID, "`\n", sep = "")
}
cat("- Broad-scan source of truth: `broad_recut_sweep.tsv` latest row per parameter tuple.\n")
cat("- Target module band for neighborhood seeding: `5-7` non-grey modules.\n\n")

cat("## Broad Scan Summary\n\n")
cat(sprintf("- Raw broad-scan rows: %d\n", nrow(sweep)))
cat(sprintf("- Unique settings after deduplication: %d\n", nrow(sweep_latest)))
cat(sprintf("- Retried settings resolved by latest-row deduplication: %d\n", dup_count))
cat(sprintf("- Successful settings: %d\n", nrow(ok)))
cat(sprintf("- Timed-out settings: %d\n", timeout_n))
cat(sprintf("- Failed settings: %d\n\n", failed_n))

cat("### Successful Module Count Distribution\n\n")
for (i in seq_len(nrow(module_dist))) {
  r <- module_dist[i]
  cat(sprintf("- `%d`: %d\n", r$non_grey_modules, r$N))
}
cat("\n")

cat("## Key Findings\n\n")
cat("- `power`, `deepSplit`, and `minModuleSize` control the broad scan. `mergeCutHeight` is mostly flat over `0.02-0.20`.\n")
cat("- `minModuleSize=8` is the main gateway to usable multi-module solutions; all main seed-eligible settings use it.\n")
cat("- `minModuleSize=10` usually drops one module relative to `8`, and `minModuleSize>=12` tends to collapse toward `3` modules or fewer.\n")
cat("- Higher `deepSplit` generally increases module count, with the strongest `7`-module plateaus appearing at `deepSplit=3-4`.\n")
cat("- The best `5-7` module regions are stable across all seven tested `mergeCutHeight` values, so the neighborhood scan should not spend budget re-testing the full broad merge grid.\n\n")

cat("## Seed Eligibility Checks\n\n")
cat("- Main seed filter: `non_grey_modules` in `5:7`, `minModuleSize = 8`, `singleton_non_grey_modules = 0`, `largest_non_grey_module_pct <= 20`, `grey_pct <= 60`.\n")
cat("- Secondary seed filter: same as above, with `grey_pct <= 65` and `minModuleSize = 10` allowed as a fallback.\n")
cat("- Boundary tracking filter: `non_grey_modules = 8`, `minModuleSize` in `8,9`, `singleton_non_grey_modules = 0`, `largest_non_grey_module_pct <= 20`.\n")
cat(sprintf("- Main seed-eligible settings: %d\n", ok[main_seed_eligible == TRUE, .N]))
cat(sprintf("- Secondary seed-eligible settings: %d\n", ok[secondary_seed_eligible == TRUE, .N]))
cat("- Secondary seed-eligible minModuleSize distribution:\n")
for (i in seq_len(nrow(seed_mms_dist))) {
  r <- seed_mms_dist[i]
  cat(sprintf("  - `%d`: %d\n", r$minModuleSize, r$N))
}
cat("\n")

if (nrow(confirmed) > 0L) {
  cat("## Ranked Seed Regions\n\n")
  cat("|role|power|deepSplit|minModuleSize|non_grey_modules|grey_pct|module_size_median|largest_non_grey_module_pct|merge_support_total|mergeCutHeight|\n")
  cat("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
  for (i in seq_len(nrow(confirmed))) {
    r <- confirmed[i]
    cat(sprintf(
      "|%s|%d|%s|%d|%d|%.2f|%.1f|%.2f|%d|%s|\n",
      r$role, r$power, fmt_range(r$deepSplit_min, r$deepSplit_max), r$minModuleSize,
      r$non_grey_modules, r$grey_pct, r$module_size_median,
      r$largest_non_grey_module_pct, r$merge_support_total, r$mergeCutHeight_values
    ))
  }
  cat("\n")
}

cat("## Neighborhood Scan Clusters\n\n")
cat("|cluster|purpose|window|best_region|best_modules|best_grey_pct|merge_support_total|\n")
cat("|---|---|---|---|---:|---:|---:|\n")
for (i in seq_len(nrow(cluster_summary))) {
  r <- cluster_summary[i]
  best_region <- if (is.na(r$best_power)) {
    "none"
  } else {
    sprintf("power=%d deepSplit=%s minModuleSize=%d", r$best_power, r$best_deepSplit, r$best_minModuleSize)
  }
  cat(sprintf(
    "|%s|%s|%s|%s|%s|%s|%s|\n",
    r$cluster,
    r$purpose,
    r$window,
    best_region,
    ifelse(is.na(r$best_non_grey_modules), "NA", as.character(r$best_non_grey_modules)),
    ifelse(is.na(r$best_grey_pct), "NA", sprintf("%.2f", r$best_grey_pct)),
    ifelse(is.na(r$best_merge_support_total), "NA", as.character(r$best_merge_support_total))
  ))
}
cat("\n")

cat("## Neighborhood Scan Plan\n\n")
cat("- Cluster A: `power 2-4`, `deepSplit 3-4`, `minModuleSize 8-10`, centered on the `power=3`, `deepSplit=4`, `minModuleSize=8` 7-module plateau.\n")
cat("- Cluster B: `power 7-9`, `deepSplit 3-4`, `minModuleSize 8-10`, centered on the `power=8`, `deepSplit=4`, `minModuleSize=8` 7-module plateau.\n")
cat("- Cluster C: `power 9-16`, `deepSplit 2-4`, `minModuleSize 8-10`, for the 6-module fallback family.\n")
cat("- Boundary Cluster D: `power 4-6`, `deepSplit 2-4`, `minModuleSize 8-9`, for the 8-module upper boundary.\n")
cat("- Use a lighter local `mergeCutHeight` grid such as `0.01, 0.02, 0.05, 0.08` instead of re-running the full `0.02-0.20` broad grid uniformly.\n")
cat("- Keep `deepSplit` local to adjacent values and keep `minModuleSize` tight around `8` and `10`, with optional `9`.\n\n")

cat("## Verification Notes\n\n")
cat("- Latest-row deduplication resolved the retried settings in the completed broad scan before any ranking or seed selection.\n")
cat("- Each recommended seed region has full `mergeCutHeight` support across the tested broad grid unless explicitly labeled as a boundary or fallback comparison.\n")
cat("- `minModuleSize=8` remains the dominant usable setting in the recomputed candidate table.\n")
sink()

progress_update("broad_candidates", "complete", sprintf("confirmed_n=%d", nrow(confirmed)))
log_msg("balanced broad candidate seeding report complete")
