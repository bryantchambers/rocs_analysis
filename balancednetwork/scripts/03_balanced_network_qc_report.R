#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA_character_) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}
run_mode <- arg_val("mode", "build")

OUT_BASE <- here("balancednetwork", "results")
OUT_TAB <- file.path(OUT_BASE, "tables")
OUT_WGCNA <- file.path(OUT_BASE, "wgcna")
OUT_STAB <- file.path(OUT_BASE, "stability")
OUT_REP <- file.path(OUT_BASE, "reports")
dir.create(OUT_REP, recursive = TRUE, showWarnings = FALSE)

sum_design <- fread(file.path(OUT_TAB, "balance_design_summary.tsv"))
bin_quota <- fread(file.path(OUT_TAB, "balance_bin_quotas.tsv"))
mods <- fread(file.path(OUT_WGCNA, "module_assignments.tsv"))
pres <- fread(file.path(OUT_WGCNA, "preservation.tsv"))
conc <- fread(file.path(OUT_WGCNA, "eigengene_concordance_age_aligned.tsv"))
stab_sum <- fread(file.path(OUT_STAB, "module_stability_summary.tsv"))
stab_run <- fread(file.path(OUT_STAB, "stability_run_summary.tsv"))

overlap_file <- file.path(OUT_STAB, "balanced_vs_original_module_overlap.tsv")
has_overlap <- file.exists(overlap_file)
if (has_overlap) overlap <- fread(overlap_file)

module_counts <- mods[, .N, by = module][order(-N)]
n_non_grey <- mods[module != "grey", uniqueN(module)]
grey_pct <- 100 * mods[module == "grey", .N] / nrow(mods)

bio_pres <- pres[module_type == "biological"]
strong_n <- bio_pres[preserved == "strong", .N]
moderate_n <- bio_pres[preserved == "moderate", .N]
weak_n <- bio_pres[preserved == "weak", .N]

mean_pearson <- mean(conc$pearson_r, na.rm = TRUE)
mean_spearman <- mean(conc$spearman_rho, na.rm = TRUE)
mean_rmse <- mean(conc$rmse, na.rm = TRUE)
mean_boot_jaccard <- mean(stab_sum$jaccard_median, na.rm = TRUE)

report <- file.path(OUT_REP, "BALANCED_NETWORK_QC_REPORT.md")
sink(report)

cat("# Balanced Network QC Report\n\n")
cat(sprintf("- Date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("- Mode: `%s`\n\n", run_mode))

cat("## 1) Balanced Design Summary\n\n")
cat("|metric|value|\n|---|---:|\n")
for (i in seq_len(nrow(sum_design))) {
  cat(sprintf("|%s|%s|\n", sum_design$metric[i], sum_design$value[i]))
}
cat("\n")

cat("Retained age-bin quotas:\n\n")
cat("|age_bin|ST8|ST13|GeoB25202_R1|quota_per_core|\n|---|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(bin_quota))) {
  cat(sprintf(
    "|%s|%d|%d|%d|%d|\n",
    bin_quota$age_bin[i], bin_quota$ST8[i], bin_quota$ST13[i], bin_quota$GeoB25202_R1[i], bin_quota$quota_per_core[i]
  ))
}
cat("\n")

cat("## 2) Balanced WGCNA Output\n\n")
cat(sprintf("- Non-grey modules: `%d`\n", n_non_grey))
cat(sprintf("- Grey fraction: `%.2f%%`\n\n", grey_pct))

cat("|module|taxa_n|\n|---|---:|\n")
for (i in seq_len(nrow(module_counts))) {
  cat(sprintf("|%s|%d|\n", module_counts$module[i], module_counts$N[i]))
}
cat("\n")

cat("## 3) Validation and Concordance\n\n")
cat(sprintf("- Biological preservation: `%d strong`, `%d moderate`, `%d weak`\n", strong_n, moderate_n, weak_n))
cat(sprintf("- Mean age-aligned Pearson: `%.3f`\n", mean_pearson))
cat(sprintf("- Mean age-aligned Spearman: `%.3f`\n", mean_spearman))
cat(sprintf("- Mean age-aligned RMSE: `%.3f`\n\n", mean_rmse))

cat("## 4) Balanced Bootstrap Stability\n\n")
cat("|module|median_jaccard|p05|p95|\n|---|---:|---:|---:|\n")
for (i in seq_len(nrow(stab_sum))) {
  cat(sprintf("|%s|%.3f|%.3f|%.3f|\n", stab_sum$module[i], stab_sum$jaccard_median[i], stab_sum$jaccard_p05[i], stab_sum$jaccard_p95[i]))
}
cat("\n")
cat(sprintf("- Mean module median Jaccard: `%.3f`\n\n", mean_boot_jaccard))

cat("|run_metric|value|\n|---|---:|\n")
for (i in seq_len(nrow(stab_run))) {
  cat(sprintf("|%s|%s|\n", stab_run$metric[i], stab_run$value[i]))
}
cat("\n")

cat("## 5) Comparison to Original exp3\n\n")
if (has_overlap) {
  cat("|balanced_module|best_original_module|best_jaccard|\n|---|---|---:|\n")
  for (i in seq_len(nrow(overlap))) {
    cat(sprintf("|%s|%s|%.3f|\n", overlap$balanced_module[i], overlap$best_original_module[i], overlap$best_jaccard[i]))
  }
  cat("\n")
} else {
  cat("Original exp3 module assignment file not found at `results/stage1/wgcna/module_assignments.tsv`.\n\n")
}

cat("## 6) Interpretation Guardrail\n\n")
cat("Use this balanced branch as a sensitivity and stress test for site/age dominance.\n")
cat("Promote into the main workflow only if stability, preservation, and biological interpretability are competitive with the current exp3 baseline.\n")

sink()

message("[balancednetwork] wrote report: ", report)
