#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

OUT_BASE <- here("networkQC", "results")
OUT_TABLE <- file.path(OUT_BASE, "tables")
OUT_REPORT <- file.path(OUT_BASE, "NETWORK_QC_DECISION_REPORT.md")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)

norm01 <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  z <- (x - rng[1]) / (rng[2] - rng[1])
  if (!higher_better) z <- 1 - z
  z
}

sweep <- fread(file.path(OUT_TABLE, "qc_parameter_sweep_summary.tsv"))
pairc <- fread(file.path(OUT_TABLE, "qc_pairwise_core_eigengene_concordance.tsv"))
setnames(pairc, old = names(pairc), new = sub("\\.V1$", "", names(pairc)))
pairc <- pairc[!is.na(pearson_r) & !is.na(spearman_rho) & !is.na(rmse)]

if (nrow(sweep[status == "ok"]) == 0) stop("No successful sweep rows found.")

# Global concordance benchmark (same for all settings currently; kept explicit for future per-setting reruns)
conc_summary <- data.table(
  pearson_mean = mean(pairc$pearson_r, na.rm = TRUE),
  spearman_mean = mean(pairc$spearman_rho, na.rm = TRUE),
  rmse_mean = mean(pairc$rmse, na.rm = TRUE)
)

dt <- sweep[status == "ok"][, .(
  power, deepSplit, mergeCutHeight, minModuleSize,
  non_grey_modules, grey_pct, module_size_median
)]

# Target assumptions for this project:
# - low grey burden
# - around 5 non-grey modules (current biologically interpreted solution)
# - avoid over-fragmentation and overly tiny modules
dt[, module_count_distance := abs(non_grey_modules - 5)]

dt[, score_grey := norm01(grey_pct, higher_better = FALSE)]
dt[, score_module_count := norm01(module_count_distance, higher_better = FALSE)]
dt[, score_module_size := norm01(module_size_median, higher_better = TRUE)]

# Static concordance score included for transparency and future extensibility.
dt[, score_concordance := mean(c(
  norm01(conc_summary$pearson_mean, TRUE),
  norm01(conc_summary$spearman_mean, TRUE),
  norm01(conc_summary$rmse_mean, FALSE)
), na.rm = TRUE)]

# Weighted composite
w_grey <- 0.45
w_count <- 0.35
w_size <- 0.20

dt[, decision_score := w_grey * score_grey +
     w_count * score_module_count +
     w_size * score_module_size]

setorder(dt, -decision_score, grey_pct, module_count_distance, -module_size_median)
dt[, rank := .I]

fwrite(dt, file.path(OUT_TABLE, "qc_decision_matrix.tsv"), sep = "\t")
fwrite(dt[1:3], file.path(OUT_TABLE, "qc_decision_top3.tsv"), sep = "\t")

sink(OUT_REPORT)
cat("# Network QC Decision Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Inputs: parameter sweep + cross-core eigengene concordance\n\n")

cat("## Scoring Strategy\n\n")
cat("The decision score is weighted toward lower grey burden and realistic module count,\n")
cat("with module-size sanity as a secondary factor.\n\n")
cat("- `score_grey` (45%): lower `grey_pct` is better\n")
cat("- `score_module_count` (35%): closer to 5 non-grey modules is better\n")
cat("- `score_module_size` (20%): larger median non-grey module size is better\n\n")

cat("## Cross-core Concordance Benchmark\n\n")
cat(sprintf("- Mean Pearson r: %.3f\n", conc_summary$pearson_mean))
cat(sprintf("- Mean Spearman rho: %.3f\n", conc_summary$spearman_mean))
cat(sprintf("- Mean RMSE: %.3f\n\n", conc_summary$rmse_mean))

cat("## Top 3 Parameter Sets\n\n")
cat("|rank|power|deepSplit|mergeCutHeight|minModuleSize|non_grey_modules|grey_pct|module_size_median|decision_score|\n")
cat("|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(min(3, nrow(dt)))) {
  r <- dt[i]
  cat(sprintf("|%d|%d|%d|%.2f|%d|%d|%.2f|%.1f|%.3f|\n",
              r$rank, r$power, r$deepSplit, r$mergeCutHeight, r$minModuleSize,
              r$non_grey_modules, r$grey_pct, r$module_size_median, r$decision_score))
}
cat("\n")

cat("## Recommendation\n\n")
cat("Run full stability diagnostics (`02b`-style bootstrap) on these top 3 settings,\n")
cat("then pick the one with best joint behavior across:\n")
cat("1) bootstrap Jaccard stability, 2) biological-module preservation,\n")
cat("3) age-aligned eigengene concordance, and 4) downstream biological consistency.\n")
sink()

message("[networkQC] decision matrix written: ", file.path(OUT_TABLE, "qc_decision_matrix.tsv"))
message("[networkQC] decision report written: ", OUT_REPORT)

