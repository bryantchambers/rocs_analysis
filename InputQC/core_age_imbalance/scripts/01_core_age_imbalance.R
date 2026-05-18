#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(here("config.R"))
set.seed(PARAMS$seed)
options(stringsAsFactors = FALSE)

BASE <- here("InputQC", "core_age_imbalance")
IN_BASE <- here("InputQC", "results")
RARE_BASE <- here("InputQC", "rarefaction_depth_qc", "results", "tables")
OUT_TABLE <- file.path(BASE, "results", "tables")
OUT_FIG <- file.path(BASE, "results", "figures")
REPORT <- file.path(BASE, "CORE_AGE_IMBALANCE_REPORT.md")
LOG <- file.path(BASE, "results", "core_age_imbalance.log")

dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) invisible(file.remove(LOG))

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

weighted_entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (!length(p)) return(NA_real_)
  -sum(p * log(p)) / log(length(p))
}

log_msg("Loading low-detection, rarefaction, and PC tables")
low <- fread(file.path(IN_BASE, "tables", "low_detection_sample_table.tsv"))
rare <- fread(file.path(RARE_BASE, "sample_rarefaction_depth_metrics.tsv"))
pc <- fread(file.path(IN_BASE, "tables", "input_pc_scores.tsv"))[
  variant == "current_taxon_centered_log",
  .(sample, PC1, PC2, PC3)
]

dt <- merge(low, rare[, .(
  sample, total_reads_from_counts, observed_taxa, shannon_vegan, vegan_chao1,
  chao1_gap, rarefaction_slope_per_10k, qc_depth_class
)], by = "sample", all.x = TRUE)
dt <- merge(dt, pc, by = "sample", all.x = TRUE, suffixes = c("", "_pc"))
dt[, site_group := fifelse(core %in% c("GeoB25202_R1", "GeoB25202_R2"), "GeoB25202", core)]

age_min <- floor(min(dt$age_kyr, na.rm = TRUE) / 10) * 10
age_max <- ceiling(max(dt$age_kyr, na.rm = TRUE) / 10) * 10
breaks_10k <- seq(age_min, age_max, by = 10)
dt[, age_bin_10k := cut(age_kyr, breaks = breaks_10k, include.lowest = TRUE, right = FALSE)]
dt[, focus_age_bin := fcase(
  age_kyr >= 40 & age_kyr <= 60, "40-60 kya",
  age_kyr >= 90 & age_kyr <= 110, "90-110 kya",
  default = "other"
)]
qbreaks <- unique(as.numeric(quantile(dt$age_kyr, probs = seq(0, 1, length.out = 7), na.rm = TRUE)))
dt[, age_bin_quantile := cut(age_kyr, breaks = qbreaks, include.lowest = TRUE)]

dt[, sequencing_attention := qc_depth_class %in% c("likely_undersequenced", "insufficient_for_call", "mixed_or_uncertain")]

log_msg("Computing age/core imbalance diagnostics")
coverage_10k <- dt[, .(
  n_samples = .N,
  median_reads = as.numeric(median(total_reads_from_counts, na.rm = TRUE)),
  median_observed_taxa = as.numeric(median(observed_taxa, na.rm = TRUE)),
  median_shannon = as.numeric(median(shannon_vegan, na.rm = TRUE)),
  median_chao1 = as.numeric(median(vegan_chao1, na.rm = TRUE)),
  low_detection_n = sum(low_reads_or_taxa %in% TRUE, na.rm = TRUE),
  low_detection_pct = mean(low_reads_or_taxa %in% TRUE, na.rm = TRUE) * 100,
  sequencing_attention_n = sum(sequencing_attention %in% TRUE, na.rm = TRUE),
  sequencing_attention_pct = mean(sequencing_attention %in% TRUE, na.rm = TRUE) * 100,
  median_PC1 = median(PC1, na.rm = TRUE),
  median_PC2 = median(PC2, na.rm = TRUE)
), by = .(age_bin_10k, core)][order(age_bin_10k, core)]

bin_totals <- dt[, .(
  bin_n = .N,
  bin_low_detection_n = sum(low_reads_or_taxa %in% TRUE, na.rm = TRUE),
  bin_attention_n = sum(sequencing_attention %in% TRUE, na.rm = TRUE)
), by = age_bin_10k]

coverage_10k <- merge(coverage_10k, bin_totals, by = "age_bin_10k", all.x = TRUE)
coverage_10k[, `:=`(
  core_share_in_bin = n_samples / bin_n,
  expected_low_detection_n = bin_low_detection_n * n_samples / bin_n,
  expected_attention_n = bin_attention_n * n_samples / bin_n
)]
coverage_10k[, `:=`(
  low_detection_enrichment = fifelse(expected_low_detection_n > 0, low_detection_n / expected_low_detection_n, NA_real_),
  attention_enrichment = fifelse(expected_attention_n > 0, sequencing_attention_n / expected_attention_n, NA_real_)
)]

coverage_10k[, fairness_label := fcase(
  bin_n < 5 | n_samples < 2, "too_sparse",
  low_detection_enrichment >= 1.5 & low_detection_pct >= 25, "enriched_beyond_sampling",
  core_share_in_bin >= 0.70 & low_detection_enrichment < 1.5, "expected_by_sampling",
  low_detection_pct > 0 & low_detection_enrichment < 1.5, "mixed_or_background",
  default = "no_low_detection_signal"
)]

focus_summary <- dt[focus_age_bin != "other", .(
  n_samples = .N,
  median_reads = as.numeric(median(total_reads_from_counts, na.rm = TRUE)),
  median_observed_taxa = as.numeric(median(observed_taxa, na.rm = TRUE)),
  median_shannon = as.numeric(median(shannon_vegan, na.rm = TRUE)),
  low_detection_n = sum(low_reads_or_taxa %in% TRUE, na.rm = TRUE),
  low_detection_pct = mean(low_reads_or_taxa %in% TRUE, na.rm = TRUE) * 100,
  sequencing_attention_n = sum(sequencing_attention %in% TRUE, na.rm = TRUE),
  sequencing_attention_pct = mean(sequencing_attention %in% TRUE, na.rm = TRUE) * 100
), by = .(focus_age_bin, core)][order(focus_age_bin, core)]

core_density <- dt[, .(
  n_samples = .N,
  age_min = min(age_kyr, na.rm = TRUE),
  age_max = max(age_kyr, na.rm = TRUE),
  age_span_kyr = max(age_kyr, na.rm = TRUE) - min(age_kyr, na.rm = TRUE),
  median_spacing_kyr = median(diff(sort(age_kyr)), na.rm = TRUE),
  samples_per_10k = .N / ((max(age_kyr, na.rm = TRUE) - min(age_kyr, na.rm = TRUE)) / 10),
  low_detection_pct = mean(low_reads_or_taxa %in% TRUE, na.rm = TRUE) * 100,
  sequencing_attention_pct = mean(sequencing_attention %in% TRUE, na.rm = TRUE) * 100
), by = core][order(core)]

bin_balance <- dt[, .N, by = .(age_bin_10k, core)]
bin_balance[, bin_n := sum(N), by = age_bin_10k]
bin_balance[, core_share := N / bin_n]
balance_summary <- bin_balance[, .(
  bin_n = first(bin_n),
  n_cores_present = .N,
  max_core_share = max(core_share),
  effective_core_balance = 1 / sum(core_share^2),
  normalized_core_entropy = weighted_entropy(core_share)
), by = age_bin_10k][order(age_bin_10k)]

winnow_steps <- data.table(
  step = c(
    "All raw-count samples",
    "Has age/core metadata",
    "Pass reads >= 50k",
    "Pass reads >= 100k",
    "Pass taxa >= 150",
    "Pass taxa >= 250",
    "Pass reads >= 50k and taxa >= 150",
    "Not low-detection",
    "Not sequencing-attention"
  ),
  n = c(
    nrow(dt),
    sum(!is.na(dt$age_kyr) & !is.na(dt$core)),
    sum(dt$total_reads_from_counts >= 50000, na.rm = TRUE),
    sum(dt$total_reads_from_counts >= 100000, na.rm = TRUE),
    sum(dt$observed_taxa >= 150, na.rm = TRUE),
    sum(dt$observed_taxa >= 250, na.rm = TRUE),
    sum(dt$total_reads_from_counts >= 50000 & dt$observed_taxa >= 150, na.rm = TRUE),
    sum(!(dt$low_reads_or_taxa %in% TRUE), na.rm = TRUE),
    sum(!(dt$sequencing_attention %in% TRUE), na.rm = TRUE)
  )
)
winnow_steps[, pct_of_total := n / first(n) * 100]
winnow_steps[, step := factor(step, levels = step)]

winnow_by_core <- rbindlist(lapply(levels(winnow_steps$step), function(step_name) {
  keep <- switch(
    as.character(step_name),
    "All raw-count samples" = rep(TRUE, nrow(dt)),
    "Has age/core metadata" = !is.na(dt$age_kyr) & !is.na(dt$core),
    "Pass reads >= 50k" = dt$total_reads_from_counts >= 50000,
    "Pass reads >= 100k" = dt$total_reads_from_counts >= 100000,
    "Pass taxa >= 150" = dt$observed_taxa >= 150,
    "Pass taxa >= 250" = dt$observed_taxa >= 250,
    "Pass reads >= 50k and taxa >= 150" = dt$total_reads_from_counts >= 50000 & dt$observed_taxa >= 150,
    "Not low-detection" = !(dt$low_reads_or_taxa %in% TRUE),
    "Not sequencing-attention" = !(dt$sequencing_attention %in% TRUE),
    rep(TRUE, nrow(dt))
  )
  dt[keep %in% TRUE, .N, by = core][, step := step_name][]
}), fill = TRUE)
winnow_by_core[, step := factor(step, levels = levels(winnow_steps$step))]

st8_rows <- coverage_10k[core == "ST8" & low_detection_n > 0]
st8_enriched_bins <- st8_rows[fairness_label == "enriched_beyond_sampling", .N]
st8_expected_bins <- st8_rows[fairness_label == "expected_by_sampling", .N]

fwrite(dt, file.path(OUT_TABLE, "core_age_sample_table.tsv"), sep = "\t")
fwrite(coverage_10k, file.path(OUT_TABLE, "age10k_core_imbalance_summary.tsv"), sep = "\t")
fwrite(focus_summary, file.path(OUT_TABLE, "focus_age_band_core_summary.tsv"), sep = "\t")
fwrite(core_density, file.path(OUT_TABLE, "core_age_density_summary.tsv"), sep = "\t")
fwrite(balance_summary, file.path(OUT_TABLE, "age10k_balance_summary.tsv"), sep = "\t")
fwrite(winnow_steps, file.path(OUT_TABLE, "sample_winnowing_summary.tsv"), sep = "\t")
fwrite(winnow_by_core, file.path(OUT_TABLE, "sample_winnowing_by_core.tsv"), sep = "\t")

log_msg("Writing imbalance figures")
core_cols <- c(GeoB25202_R1 = "#1b9e77", GeoB25202_R2 = "#66a61e", ST8 = "#7570b3", ST13 = "#d95f02")

p_waterfall <- ggplot(winnow_steps, aes(step, n)) +
  geom_col(fill = "#4c78a8", width = 0.72) +
  geom_text(aes(label = sprintf("%d\n%.0f%%", n, pct_of_total)), vjust = -0.25, size = 3) +
  labs(title = "Sample winnowing by candidate QC rules", x = NULL, y = "Samples retained") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(OUT_FIG, "sample_winnowing_waterfall.png"), p_waterfall, width = 10.5, height = 5.5, dpi = 160)

p_winnow_core <- ggplot(winnow_by_core, aes(step, N, fill = core)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = core_cols) +
  labs(title = "Sample winnowing by core", x = NULL, y = "Samples retained", fill = "Core") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(OUT_FIG, "sample_winnowing_by_core.png"), p_winnow_core, width = 10.5, height = 5.5, dpi = 160)

p_coverage <- ggplot(coverage_10k, aes(age_bin_10k, core, fill = n_samples)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_samples), size = 3) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  labs(title = "Core coverage by 10 kyr age bin", x = "Age bin", y = NULL, fill = "n") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "core_age_coverage_heatmap.png"), p_coverage, width = 11, height = 4.8, dpi = 160)

p_low <- ggplot(coverage_10k, aes(age_bin_10k, core, fill = low_detection_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", low_detection_pct)), size = 3) +
  scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey90") +
  labs(title = "Low-detection fraction by core and age bin", x = "Age bin", y = NULL, fill = "Low-detection") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "low_detection_fraction_heatmap.png"), p_low, width = 11, height = 4.8, dpi = 160)

p_attn <- ggplot(coverage_10k, aes(age_bin_10k, core, fill = sequencing_attention_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", sequencing_attention_pct)), size = 3) +
  scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey90") +
  labs(title = "Sequencing-attention fraction by core and age bin", x = "Age bin", y = NULL, fill = "Attention") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "sequencing_attention_fraction_heatmap.png"), p_attn, width = 11, height = 4.8, dpi = 160)

p_enrich <- ggplot(coverage_10k[!is.na(low_detection_enrichment)], aes(age_bin_10k, core, fill = low_detection_enrichment)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", low_detection_enrichment)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 1, na.value = "grey90") +
  labs(title = "Low-detection enrichment vs age-bin expectation", x = "Age bin", y = NULL, fill = "Observed / expected") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "low_detection_enrichment_heatmap.png"), p_enrich, width = 11, height = 4.8, dpi = 160)

p_density <- ggplot(core_density, aes(core, samples_per_10k, fill = core)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f", samples_per_10k)), vjust = -0.25, size = 3) +
  scale_fill_manual(values = core_cols) +
  labs(title = "Sampling density by core", x = NULL, y = "Samples per 10 kyr") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

infographic <- (p_waterfall / p_winnow_core / p_density) + plot_layout(heights = c(1.1, 1.1, 0.8))
ggsave(file.path(OUT_FIG, "sample_winnowing_infographic.png"), infographic, width = 11, height = 13, dpi = 160)

sink(REPORT)
cat("# Core/Age Imbalance Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Progress log: `InputQC/core_age_imbalance/results/core_age_imbalance.log`\n\n")

cat("## Main interpretation\n\n")
cat("This analysis checks whether ST8 low-detection clustering is stronger than expected after accounting for age-bin sample availability. A cluster is more concerning when ST8 has more low-detection samples than expected from its share of samples in the same age bin.\n\n")
cat(sprintf("ST8 has `%d` age bins with low-detection samples enriched beyond sampling expectation and `%d` bins where the pattern is mostly expected from sampling dominance.\n\n",
            st8_enriched_bins, st8_expected_bins))

cat("## Core sampling density\n\n")
cat("|core|n|age min|age max|age span kyr|samples per 10 kyr|low-detection pct|sequencing-attention pct|\n")
cat("|---|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(core_density))) {
  cat(sprintf("|%s|%d|%.2f|%.2f|%.2f|%.2f|%.1f|%.1f|\n",
              core_density$core[i], core_density$n_samples[i], core_density$age_min[i],
              core_density$age_max[i], core_density$age_span_kyr[i],
              core_density$samples_per_10k[i], core_density$low_detection_pct[i],
              core_density$sequencing_attention_pct[i]))
}
cat("\n")

cat("## Focus age bands\n\n")
if (nrow(focus_summary) == 0) {
  cat("No samples fell in the focused 40-60 kya or 90-110 kya windows.\n\n")
} else {
  cat("|age band|core|n|median reads|median observed taxa|low-detection pct|sequencing-attention pct|\n")
  cat("|---|---|---:|---:|---:|---:|---:|\n")
  for (i in seq_len(nrow(focus_summary))) {
    cat(sprintf("|%s|%s|%d|%.0f|%.0f|%.1f|%.1f|\n",
                focus_summary$focus_age_bin[i], focus_summary$core[i], focus_summary$n_samples[i],
                focus_summary$median_reads[i], focus_summary$median_observed_taxa[i],
                focus_summary$low_detection_pct[i], focus_summary$sequencing_attention_pct[i]))
  }
  cat("\n")
}

cat("## How to read the fairness labels\n\n")
cat("- `expected_by_sampling`: a core dominates the age bin, so some clustering from that core is expected.\n")
cat("- `enriched_beyond_sampling`: low-detection samples from that core exceed what its age-bin sample share predicts.\n")
cat("- `mixed_or_background`: low-detection exists but is not strongly enriched beyond sampling.\n")
cat("- `too_sparse`: too few samples to make a fair comparison.\n\n")

cat("## Outputs\n\n")
cat("- Sample table: `InputQC/core_age_imbalance/results/tables/core_age_sample_table.tsv`\n")
cat("- 10 kyr imbalance summary: `InputQC/core_age_imbalance/results/tables/age10k_core_imbalance_summary.tsv`\n")
cat("- Focus age bands: `InputQC/core_age_imbalance/results/tables/focus_age_band_core_summary.tsv`\n")
cat("- Core density: `InputQC/core_age_imbalance/results/tables/core_age_density_summary.tsv`\n")
cat("- Winnowing summary: `InputQC/core_age_imbalance/results/tables/sample_winnowing_summary.tsv`\n")
cat("- Figures: `InputQC/core_age_imbalance/results/figures/`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
