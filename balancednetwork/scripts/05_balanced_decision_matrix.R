#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(here::here("balancednetwork", "config_balanced.R"))

sweep <- fread(file.path(BAL$qc_tables_dir, "qc_parameter_sweep_summary.tsv"))
ok <- sweep[status == "ok"]
if (nrow(ok) == 0) stop("No successful balanced sweep rows.")

dt <- ok[, .(
  power, deepSplit, mergeCutHeight, minModuleSize,
  non_grey_modules, grey_pct, module_size_median
)]

dt[, module_count_distance := abs(non_grey_modules - 5)]
dt[, score_grey := norm01(grey_pct, higher_better = FALSE)]
dt[, score_module_count := norm01(module_count_distance, higher_better = FALSE)]
dt[, score_module_size := norm01(module_size_median, higher_better = TRUE)]

w_grey <- 0.45
w_count <- 0.35
w_size <- 0.20

dt[, decision_score := w_grey * score_grey + w_count * score_module_count + w_size * score_module_size]
setorder(dt, -decision_score, grey_pct, module_count_distance, -module_size_median)
dt[, rank := .I]

fwrite(dt, file.path(BAL$qc_tables_dir, "qc_decision_matrix.tsv"), sep = "\t")
fwrite(dt[1:min(5, .N)], file.path(BAL$qc_tables_dir, "qc_decision_top5.tsv"), sep = "\t")

report_file <- file.path(BAL$qc_dir, "BALANCED_QC_DECISION_REPORT.md")
sink(report_file)
cat("# Balanced Network QC Decision Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
cat("## Scoring Strategy\n\n")
cat("- `score_grey` (45%): lower `grey_pct` is better\n")
cat("- `score_module_count` (35%): closer to 5 non-grey modules is better\n")
cat("- `score_module_size` (20%): larger median non-grey module size is better\n\n")
cat("## Top 5 Parameter Sets\n\n")
cat("|rank|power|deepSplit|mergeCutHeight|minModuleSize|non_grey_modules|grey_pct|module_size_median|decision_score|\n")
cat("|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
top5 <- dt[1:min(5, .N)]
for (i in seq_len(nrow(top5))) {
  r <- top5[i]
  cat(sprintf("|%d|%d|%d|%.2f|%d|%d|%.2f|%.1f|%.3f|\n",
              r$rank, r$power, r$deepSplit, r$mergeCutHeight, r$minModuleSize,
              r$non_grey_modules, r$grey_pct, r$module_size_median, r$decision_score))
}
cat("\n")
cat("## Next Step\n\n")
cat("Run full balanced evaluation on the top 5 settings with high bootstrap depth and rank by multi-metric final score.\n")
sink()

log_msg("balanced decision matrix complete")
