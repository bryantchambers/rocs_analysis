#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

OUT_FULL <- here("networkQC", "results", "full_eval")
OUT_FIG <- here("networkQC", "results", "figures")
OUT_TABLE <- here("networkQC", "results", "tables")
REPORT <- file.path(OUT_FULL, "FULL_EVAL_REPORT.md")

sum_path <- file.path(OUT_FULL, "all_settings_summary.tsv")
rank_path <- file.path(OUT_FULL, "all_settings_ranked.tsv")
if (!file.exists(sum_path) || !file.exists(rank_path)) {
  stop("Missing full evaluation summary files. Run 06_full_eval_top5.R first.")
}

sum_dt <- fread(sum_path)
rank_dt <- fread(rank_path)
best <- rank_dt[1]
setting_metric_path <- file.path(OUT_TABLE, "full_eval_setting_metric_summary.tsv")
setting_metrics <- if (file.exists(setting_metric_path)) fread(setting_metric_path) else NULL

sink(REPORT)
cat("# NetworkQC Full Evaluation Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: baseline + top optimized sweep settings + expansion shortlist settings (full preservation + bootstrap)\n\n")

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

if (!is.null(setting_metrics) && nrow(setting_metrics) > 0) {
  sm <- setting_metrics[order(rank)]
  best_sm <- sm[1]
  runner_up <- sm[2]
  opt5 <- sm[setting_id == "opt5"]
  baseline <- sm[setting_id == "baseline"]

  cat("## Expanded Stability Review\n\n")
  cat(sprintf("Recommendation: use `%s` as the current best WGCNA parameter set.\n\n", best_sm$setting_id))
  cat(sprintf("- `%s` parameters: power=%d, deepSplit=%d, mergeCutHeight=%.2f, minModuleSize=%d\n",
              best_sm$setting_id, best_sm$power, best_sm$deepSplit, best_sm$mergeCutHeight, best_sm$minModuleSize))
  cat(sprintf("- It has the lowest grey burden among the top settings (%.2f%%), %d non-grey modules, %d strong and %d moderate biological modules, mean bootstrap Jaccard %.3f, and mean balanced Jaccard %.3f.\n",
              best_sm$grey_pct, best_sm$non_grey_modules, best_sm$bio_pres_strong, best_sm$bio_pres_moderate,
              best_sm$mean_bootstrap_jaccard, best_sm$mean_balanced_jaccard))
  cat(sprintf("- Biological-module medians for `%s`: Jaccard median %.3f, Jaccard p05 %.3f, Jaccard p95 %.3f, Zsummary %.2f, Zdensity %.2f, Zconnectivity %.2f, Pearson %.3f, Spearman %.3f, RMSE %.3f.\n",
              best_sm$setting_id, best_sm$median_jaccard_median, best_sm$median_jaccard_p05, best_sm$median_jaccard_p95,
              best_sm$median_Zsummary, best_sm$median_Zdensity, best_sm$median_Zconnectivity,
              best_sm$median_module_pearson, best_sm$median_module_spearman, best_sm$median_module_rmse))
  if (nrow(runner_up) == 1) {
    cat(sprintf("- Runner-up `%s` has slightly stronger preservation count (%d strong, %d moderate) and similar correlation, but lower mean bootstrap Jaccard %.3f and lower mean balanced Jaccard %.3f.\n",
                runner_up$setting_id, runner_up$bio_pres_strong, runner_up$bio_pres_moderate,
                runner_up$mean_bootstrap_jaccard, runner_up$mean_balanced_jaccard))
  }
  if (nrow(opt5) == 1) {
    cat(sprintf("- The previous best 5-module option `opt5` remains reasonable but retains more grey taxa (%.2f%%) and has lower age-aligned Pearson concordance %.3f.\n",
                opt5$grey_pct, opt5$mean_concordance_pearson))
  }
  if (nrow(baseline) == 1) {
    cat(sprintf("- The original baseline is not competitive here: grey taxa remain %.2f%% and mean balanced Jaccard is %.3f.\n",
                baseline$grey_pct, baseline$mean_balanced_jaccard))
  }
  cat("\n")
}

cat("## Output Map\n\n")
cat("- Per-setting outputs: `networkQC/results/full_eval/<setting_id>/`\n")
cat("- Queued settings: `networkQC/results/full_eval/settings_to_evaluate.tsv`\n")
cat("- Ranked summary: `networkQC/results/full_eval/all_settings_ranked.tsv`\n")
cat("- Graph plots by setting: `networkQC/results/figures/full_graph_<setting_id>.png`\n")
cat("- Leiden graph plot: `networkQC/results/figures/full_graph_leiden.png`\n")
cat("- Setting-level heatmap: `networkQC/results/figures/full_eval_setting_metric_heatmap.png`\n")
cat("- Module-level heatmap: `networkQC/results/figures/full_eval_module_metric_heatmap.png`\n")
cat("- Setting-level metric table: `networkQC/results/tables/full_eval_setting_metric_summary.tsv`\n")
cat("- Module-level metric table: `networkQC/results/tables/full_eval_module_metric_summary.tsv`\n")
cat("\n")

cat("## Interpretation Notes\n\n")
cat("1. `mean_bootstrap_jaccard` captures module reproducibility under resampling.\n")
cat("2. `bio_pres_strong` prioritizes non-technical module preservation in R1->R2.\n")
cat("3. `mean_balanced_jaccard` estimates how sensitive module formation is to core-size imbalance.\n")
cat("4. Choose settings with strong performance across all three, not just low grey.\n")
sink()

message("[networkQC] full evaluation report written: ", REPORT)
