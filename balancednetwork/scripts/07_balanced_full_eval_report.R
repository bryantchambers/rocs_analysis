#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(here::here("balancednetwork", "config_balanced.R"))

ranked_path <- file.path(BAL$qc_full_dir, "all_settings_ranked.tsv")
if (!file.exists(ranked_path)) stop("Missing all_settings_ranked.tsv. Run full eval first.")
ranked <- fread(ranked_path)

report <- file.path(BAL$qc_dir, "BALANCED_FULL_EVAL_REPORT.md")
sink(report)

cat("# Balanced Network Full Evaluation Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
cat("## Ranking Table\n\n")
cat("|rank|setting_id|power|deepSplit|mergeCutHeight|minModuleSize|grey_pct|mean_bootstrap_jaccard|bio_pres_strong|mean_concordance_pearson|final_score|n_boot_successful|n_boot_failed|\n")
cat("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(ranked))) {
  r <- ranked[i]
  cat(sprintf("|%d|%s|%d|%d|%.2f|%d|%.2f|%.3f|%d|%.3f|%.3f|%d|%d|\n",
              r$rank, r$setting_id, r$power, r$deepSplit, r$mergeCutHeight, r$minModuleSize,
              r$grey_pct, r$mean_bootstrap_jaccard, r$bio_pres_strong, r$mean_concordance_pearson,
              r$final_score, r$n_boot_successful, r$n_boot_failed))
}
cat("\n")

best <- ranked[1]
cat("## Best Setting Recommendation\n\n")
cat(sprintf("- Best setting: `%s`\n", best$setting_id))
cat(sprintf("- Parameters: power=%d, deepSplit=%d, mergeCutHeight=%.2f, minModuleSize=%d\n",
            best$power, best$deepSplit, best$mergeCutHeight, best$minModuleSize))
cat(sprintf("- Grey fraction: %.2f%%\n", best$grey_pct))
cat(sprintf("- Mean bootstrap Jaccard: %.3f\n", best$mean_bootstrap_jaccard))
cat(sprintf("- Biological strong preservation modules: %d\n", best$bio_pres_strong))
cat(sprintf("- Mean concordance Pearson: %.3f\n", best$mean_concordance_pearson))
cat(sprintf("- Successful bootstraps: %d (failed: %d)\n", best$n_boot_successful, best$n_boot_failed))
cat("\n")

cat("## Notes\n\n")
cat("This ranking is balanced-network-specific and should be compared against the original unbalanced exp3 branch before promotion.\n")
sink()

log_msg("balanced full eval report complete")
