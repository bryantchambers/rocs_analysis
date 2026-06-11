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
sweep_file <- file.path(BAL$qc_broad_dir, "broad_recut_sweep.tsv")
if (!file.exists(sweep_file)) stop("Missing broad recut sweep: ", sweep_file)

progress_update("broad_candidates", "start")

sweep <- fread(sweep_file)
ok <- sweep[status == "ok"]
if (nrow(ok) == 0) stop("No successful broad recut rows.")

ok[, setting_id := sprintf("broad_%03d", .I)]
ok[, module_target_distance := abs(non_grey_modules - BAL_PARAMS$broad_target_modules)]
ok[, keep_candidate := (
  non_grey_modules >= BAL_PARAMS$broad_keep_module_min &
    non_grey_modules <= BAL_PARAMS$broad_keep_module_max &
    grey_pct <= BAL_PARAMS$broad_keep_grey_pct_max &
    module_size_median >= BAL_PARAMS$broad_keep_module_median_min &
    largest_non_grey_module_pct <= BAL_PARAMS$broad_keep_largest_module_pct_max &
    singleton_non_grey_modules == 0
)]

candidates <- ok[keep_candidate == TRUE][
  order(module_target_distance, grey_pct, -module_size_median, largest_non_grey_module_pct)
]

fwrite(candidates, file.path(BAL$qc_broad_dir, "broad_recut_candidates.tsv"), sep = "\t")
progress_update("broad_candidates", "ok", sprintf("n_candidates=%d", nrow(candidates)))

selected <- candidates[1:min(BAL_PARAMS$broad_top_n_confirm, .N)]
confirm_rows <- vector("list", nrow(selected))

for (i in seq_len(nrow(selected))) {
  pars <- selected[i]
  sid <- pars$setting_id
  progress_update(
    "broad_confirm_setting",
    "start",
    sprintf(
      "%s power=%d deepSplit=%d mergeCutHeight=%.2f minModuleSize=%d",
      sid, pars$power, pars$deepSplit, pars$mergeCutHeight, pars$minModuleSize
    )
  )
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

  if (inherits(fit, "error")) {
    confirm_rows[[i]] <- data.table(
      setting_id = sid,
      power = pars$power,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      minModuleSize = pars$minModuleSize,
      confirm_status = "failed",
      confirm_error = conditionMessage(fit),
      confirmed_match = FALSE
    )
    progress_update("broad_confirm_setting", "failed", sprintf("%s %s", sid, conditionMessage(fit)))
    next
  }

  direct_sum <- module_summary_from_colors(fit$colors)
  confirm_rows[[i]] <- cbind(
    data.table(
      setting_id = sid,
      power = pars$power,
      deepSplit = pars$deepSplit,
      mergeCutHeight = pars$mergeCutHeight,
      minModuleSize = pars$minModuleSize,
      confirm_status = "ok",
      confirm_error = ""
    ),
    direct_sum
  )[
    ,
    confirmed_match := (
      non_grey_modules == pars$non_grey_modules &
        isTRUE(all.equal(grey_pct, pars$grey_pct, tolerance = 1e-8)) &
        isTRUE(all.equal(module_size_median, pars$module_size_median, tolerance = 1e-8)) &
        isTRUE(all.equal(largest_non_grey_module_pct, pars$largest_non_grey_module_pct, tolerance = 1e-8))
    )
  ]
  progress_update(
    "broad_confirm_setting",
    "ok",
    sprintf("%s confirmed_match=%s", sid, confirm_rows[[i]]$confirmed_match[[1]])
  )
}

confirmed <- if (length(confirm_rows)) rbindlist(confirm_rows, fill = TRUE) else data.table()
fwrite(confirmed, file.path(BAL$qc_broad_dir, "broad_recut_confirmed.tsv"), sep = "\t")

report_file <- file.path(BAL$qc_broad_dir, "BALANCED_BROAD_SWEEP_REPORT.md")
sink(report_file)
cat("# Balanced Broad Sweep Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Target module band: 5-7 non-grey modules, with acceptance widened to 4-9 for discovery.\n\n")
cat("## Candidate Filter\n\n")
cat(sprintf("- `non_grey_modules`: %d-%d\n", BAL_PARAMS$broad_keep_module_min, BAL_PARAMS$broad_keep_module_max))
cat(sprintf("- `grey_pct <= %.1f`\n", BAL_PARAMS$broad_keep_grey_pct_max))
cat(sprintf("- `module_size_median >= %.1f`\n", BAL_PARAMS$broad_keep_module_median_min))
cat(sprintf("- `largest_non_grey_module_pct <= %.1f`\n", BAL_PARAMS$broad_keep_largest_module_pct_max))
cat("- `singleton_non_grey_modules = 0`\n\n")
cat("## Broad Recut Summary\n\n")
cat(sprintf("- Successful settings: %d\n", nrow(ok)))
cat(sprintf("- Candidate settings: %d\n", nrow(candidates)))
cat(sprintf("- Confirmed settings: %d\n\n", confirmed[confirm_status == "ok" & confirmed_match == TRUE, .N]))

if (nrow(confirmed) > 0) {
  show_dt <- confirmed[confirm_status == "ok"][order(abs(non_grey_modules - BAL_PARAMS$broad_target_modules), grey_pct)][1:min(10, .N)]
  cat("## Top Confirmed Settings\n\n")
  cat("|setting_id|power|deepSplit|mergeCutHeight|minModuleSize|non_grey_modules|grey_pct|module_size_median|largest_non_grey_module_pct|confirmed_match|\n")
  cat("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
  for (i in seq_len(nrow(show_dt))) {
    r <- show_dt[i]
    cat(sprintf(
      "|%s|%d|%d|%.2f|%d|%d|%.2f|%.1f|%.2f|%s|\n",
      r$setting_id, r$power, r$deepSplit, r$mergeCutHeight, r$minModuleSize,
      r$non_grey_modules, r$grey_pct, r$module_size_median,
      r$largest_non_grey_module_pct, ifelse(isTRUE(r$confirmed_match), "yes", "no")
    ))
  }
}
sink()

progress_update("broad_candidates", "complete", sprintf("confirmed_n=%d", nrow(confirmed)))
log_msg("balanced broad candidate confirmation complete")
