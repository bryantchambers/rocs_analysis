#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

OUT_FULL <- here("networkQC", "results", "full_eval")
OUT_FIG <- here("networkQC", "results", "figures")
REPORT <- file.path(OUT_FULL, "FULL_EVAL_REPORT.md")

sum_path <- file.path(OUT_FULL, "all_settings_summary.tsv")
rank_path <- file.path(OUT_FULL, "all_settings_ranked.tsv")
if (!file.exists(sum_path) || !file.exists(rank_path)) {
  stop("Missing full evaluation summary files. Run 06_full_eval_top5.R first.")
}

sum_dt <- fread(sum_path)
rank_dt <- fread(rank_path)
best <- rank_dt[1]

sink(REPORT)
cat("# NetworkQC Full Evaluation Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: baseline vs top optimized sweep settings (full preservation + bootstrap)\n\n")

cat("## Ranking Table\n\n")
cat("|rank|setting_id|power|deepSplit|mergeCutHeight|minModuleSize|grey_pct|mean_bootstrap_jaccard|bio_pres_strong|mean_concordance_pearson|mean_balanced_jaccard|final_score|\n")
cat("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(rank_dt))) {
  r <- rank_dt[i]
  cat(sprintf("|%d|%s|%d|%d|%.2f|%d|%.2f|%.3f|%d|%.3f|%.3f|%.3f|\n",
              r$rank, r$setting_id, r$power, r$deepSplit, r$mergeCutHeight, r$minModuleSize,
              r$grey_pct, r$mean_bootstrap_jaccard, r$bio_pres_strong,
              r$mean_concordance_pearson, r$mean_balanced_jaccard, r$final_score))
}
cat("\n")

cat("## Best Setting Recommendation\n\n")
cat(sprintf("- Best setting: `%s`\n", best$setting_id))
cat(sprintf("- Parameters: power=%d, deepSplit=%d, mergeCutHeight=%.2f, minModuleSize=%d\n",
            best$power, best$deepSplit, best$mergeCutHeight, best$minModuleSize))
cat(sprintf("- Grey fraction: %.2f%%\n", best$grey_pct))
cat(sprintf("- Bootstrap stability (mean module median Jaccard): %.3f\n", best$mean_bootstrap_jaccard))
cat(sprintf("- Biological preservation strong modules: %d\n", best$bio_pres_strong))
cat(sprintf("- Core-balance robustness (mean balanced Jaccard): %.3f\n", best$mean_balanced_jaccard))
cat("\n")

cat("## Output Map\n\n")
cat("- Per-setting outputs: `networkQC/results/full_eval/<setting_id>/`\n")
cat("- Ranked summary: `networkQC/results/full_eval/all_settings_ranked.tsv`\n")
cat("- Graph plots by setting: `networkQC/results/figures/full_graph_<setting_id>.png`\n")
cat("- Leiden graph plot: `networkQC/results/figures/full_graph_leiden.png`\n")
cat("\n")

cat("## Interpretation Notes\n\n")
cat("1. `mean_bootstrap_jaccard` captures module reproducibility under resampling.\n")
cat("2. `bio_pres_strong` prioritizes non-technical module preservation in R1->R2.\n")
cat("3. `mean_balanced_jaccard` estimates how sensitive module formation is to core-size imbalance.\n")
cat("4. Choose settings with strong performance across all three, not just low grey.\n")
sink()

message("[networkQC] full evaluation report written: ", REPORT)

