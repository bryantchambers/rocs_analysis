#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(phyloseq)
  library(vegan)
  library(iNEXT)
  library(patchwork)
})

source(here("config.R"))
set.seed(PARAMS$seed)
options(stringsAsFactors = FALSE)

BASE <- here("InputQC", "rarefaction_depth_qc")
IN_BASE <- here("InputQC", "results")
OUT_TABLE <- file.path(BASE, "results", "tables")
OUT_FIG <- file.path(BASE, "results", "figures")
REPORT <- file.path(BASE, "RAREFACTION_DEPTH_QC_REPORT.md")
LOG <- file.path(BASE, "results", "rarefaction_depth_qc.log")

dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) invisible(file.remove(LOG))

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

safe_diversity <- function(x) {
  x <- as.numeric(x)
  if (!sum(x, na.rm = TRUE)) return(NA_real_)
  vegan::diversity(x, index = "shannon")
}

safe_rarefy <- function(comm, depth) {
  out <- rep(NA_real_, nrow(comm))
  ok <- rowSums(comm) >= depth & depth > 0
  if (any(ok)) {
    out[ok] <- suppressWarnings(as.numeric(vegan::rarefy(comm[ok, , drop = FALSE], sample = depth)))
  }
  out
}

safe_iNEXT_info <- function(counts) {
  abund_list <- lapply(seq_len(ncol(counts)), function(i) {
    x <- as.numeric(counts[, i])
    x[x > 0]
  })
  names(abund_list) <- colnames(counts)
  info <- tryCatch(
    as.data.table(iNEXT::DataInfo(abund_list, datatype = "abundance")),
    error = function(e) data.table(Assemblage = names(abund_list), SC = NA_real_)
  )
  setnames(info, "Assemblage", "sample", skip_absent = TRUE)
  info
}

log_msg("Loading raw counts, covariates, low-detection flags, and PC scores")
raw_counts <- readRDS(file.path(IN_BASE, "inputs", "raw_counts.rds"))
cov <- fread(file.path(IN_BASE, "tables", "sample_technical_covariates.tsv"))
low <- fread(file.path(IN_BASE, "tables", "low_detection_sample_table.tsv"))
pc <- fread(file.path(IN_BASE, "tables", "input_pc_scores.tsv"))[
  variant == "current_taxon_centered_log",
  .(sample, PC1, PC2, PC3)
]

common <- Reduce(intersect, list(colnames(raw_counts), cov$sample, low$sample))
raw_counts <- raw_counts[, common, drop = FALSE]
cov <- cov[match(common, sample)]
low <- low[match(common, sample)]
pc <- pc[match(common, sample)]

sample_df <- as.data.frame(cov)
rownames(sample_df) <- sample_df$sample
ps <- phyloseq(
  otu_table(raw_counts, taxa_are_rows = TRUE),
  sample_data(sample_df)
)
saveRDS(ps, file.path(OUT_TABLE, "phyloseq_input_object.rds"))

comm <- t(raw_counts)
log_msg("Computing diversity, Chao1, rarefaction, and coverage metrics")
phy_rich <- as.data.table(phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Chao1")))
phy_rich[, sample := rownames(phyloseq::estimate_richness(ps, measures = c("Observed")))]

est <- as.data.table(t(vegan::estimateR(comm)), keep.rownames = "sample")
setnames(est, c("S.obs", "S.chao1", "se.chao1"), c("vegan_observed", "vegan_chao1", "vegan_chao1_se"), skip_absent = TRUE)

metrics <- data.table(
  sample = rownames(comm),
  total_reads_from_counts = as.numeric(rowSums(comm)),
  observed_taxa = as.integer(rowSums(comm > 0)),
  singleton_taxa = as.integer(rowSums(comm == 1)),
  doubleton_taxa = as.integer(rowSums(comm == 2)),
  shannon_vegan = apply(comm, 1, safe_diversity)
)
metrics <- merge(metrics, phy_rich[, .(sample, phyloseq_observed = Observed, phyloseq_shannon = Shannon, phyloseq_chao1 = Chao1)], by = "sample", all.x = TRUE)
metrics <- merge(metrics, est[, .(sample, vegan_observed, vegan_chao1, vegan_chao1_se)], by = "sample", all.x = TRUE)

iinfo <- safe_iNEXT_info(raw_counts)
if ("SC" %in% names(iinfo)) {
  metrics <- merge(metrics, iinfo[, .(sample, iNEXT_sample_coverage = SC)], by = "sample", all.x = TRUE)
} else {
  metrics[, iNEXT_sample_coverage := NA_real_]
}

checkpoints <- c(5000, 10000, 25000, 50000, 100000, 250000, 500000)
rare_dt <- rbindlist(lapply(checkpoints, function(depth) {
  data.table(sample = rownames(comm), rarefaction_depth = depth, rarefied_observed = safe_rarefy(comm, depth))
}), fill = TRUE)

depth80 <- pmax(1, floor(metrics$total_reads_from_counts * 0.80))
rich80 <- vapply(seq_len(nrow(comm)), function(i) {
  if (metrics$total_reads_from_counts[i] < 10) return(NA_real_)
  suppressWarnings(as.numeric(vegan::rarefy(comm[i, , drop = FALSE], sample = depth80[i])))
}, numeric(1))

metrics[, richness_at_80pct_depth := rich80]
metrics[, `:=`(
  rarefaction_slope_per_10k = fifelse(total_reads_from_counts > depth80,
                                      (observed_taxa - richness_at_80pct_depth) / (total_reads_from_counts - depth80) * 10000,
                                      NA_real_),
  chao1_gap = fifelse(is.finite(vegan_chao1) & vegan_chao1 > 0,
                      pmax(0, (vegan_chao1 - observed_taxa) / vegan_chao1),
                      NA_real_),
  singleton_fraction = fifelse(observed_taxa > 0, singleton_taxa / observed_taxa, NA_real_),
  doubleton_fraction = fifelse(observed_taxa > 0, doubleton_taxa / observed_taxa, NA_real_)
)]

metrics <- merge(metrics, cov, by = "sample", all.x = TRUE, suffixes = c("", "_cov"))
metrics <- merge(metrics, low[, .(
  sample, low_reads_50k, low_reads_100k, low_taxa_150, low_taxa_250,
  low_reads_or_taxa, low_reads_and_taxa
)], by = "sample", all.x = TRUE)
metrics <- merge(metrics, pc, by = "sample", all.x = TRUE)
metrics[, site_group := fifelse(core %in% c("GeoB25202_R1", "GeoB25202_R2"), "GeoB25202", core)]
metrics[, focus_age_band := fifelse(age_kyr >= 40 & age_kyr <= 60, "40-60 kya",
                             fifelse(age_kyr >= 90 & age_kyr <= 110, "90-110 kya", "other"))]

# These thresholds are deliberately conservative and interpretable. They are QC flags, not automatic exclusions.
metrics[, qc_depth_class := fcase(
  total_reads_from_counts < 5000 | observed_taxa < 25, "insufficient_for_call",
  is.finite(chao1_gap) & chao1_gap >= 0.25 & is.finite(rarefaction_slope_per_10k) & rarefaction_slope_per_10k >= 10, "likely_undersequenced",
  observed_taxa < 250 & is.finite(rarefaction_slope_per_10k) & rarefaction_slope_per_10k < 5, "low_diversity_but_saturated",
  is.finite(chao1_gap) & chao1_gap < 0.10 & is.finite(rarefaction_slope_per_10k) & rarefaction_slope_per_10k < 5, "adequate_depth",
  default = "mixed_or_uncertain"
)]

metrics[, chao1_caution := vegan_chao1 == observed_taxa | singleton_taxa == 0]

threshold_summary <- rbindlist(lapply(c("low_reads_50k", "low_reads_100k", "low_taxa_150", "low_taxa_250", "low_reads_or_taxa", "low_reads_and_taxa"), function(flag) {
  metrics[, .(
    n = .N,
    pct = .N / nrow(metrics) * 100,
    median_reads = as.numeric(median(total_reads_from_counts, na.rm = TRUE)),
    median_observed_taxa = as.numeric(median(observed_taxa, na.rm = TRUE)),
    median_shannon = as.numeric(median(shannon_vegan, na.rm = TRUE)),
    median_chao1_gap = as.numeric(median(chao1_gap, na.rm = TRUE)),
    median_slope_per_10k = as.numeric(median(rarefaction_slope_per_10k, na.rm = TRUE)),
    pct_likely_undersequenced = mean(qc_depth_class == "likely_undersequenced", na.rm = TRUE) * 100,
    pct_low_diversity_saturated = mean(qc_depth_class == "low_diversity_but_saturated", na.rm = TRUE) * 100
  ), by = .(flag_value = get(flag))][, threshold_flag := flag][]
}), fill = TRUE)
setcolorder(threshold_summary, c("threshold_flag", "flag_value"))

class_summary <- metrics[, .(
  n = .N,
  pct = .N / nrow(metrics) * 100,
  median_reads = as.numeric(median(total_reads_from_counts, na.rm = TRUE)),
  median_observed_taxa = as.numeric(median(observed_taxa, na.rm = TRUE)),
  median_shannon = as.numeric(median(shannon_vegan, na.rm = TRUE)),
  median_chao1_gap = as.numeric(median(chao1_gap, na.rm = TRUE)),
  median_slope_per_10k = as.numeric(median(rarefaction_slope_per_10k, na.rm = TRUE))
), by = qc_depth_class][order(-n)]

core_class_summary <- metrics[, .N, by = .(core, qc_depth_class)][
  , pct_within_core := N / sum(N) * 100, by = core
][order(core, -N)]

fwrite(metrics, file.path(OUT_TABLE, "sample_rarefaction_depth_metrics.tsv"), sep = "\t")
fwrite(rare_dt, file.path(OUT_TABLE, "rarefaction_checkpoint_table.tsv"), sep = "\t")
fwrite(threshold_summary, file.path(OUT_TABLE, "threshold_depth_summary.tsv"), sep = "\t")
fwrite(class_summary, file.path(OUT_TABLE, "qc_depth_class_summary.tsv"), sep = "\t")
fwrite(core_class_summary, file.path(OUT_TABLE, "core_qc_depth_class_summary.tsv"), sep = "\t")

log_msg("Writing rarefaction and depth figures")
core_cols <- c(GeoB25202_R1 = "#1b9e77", GeoB25202_R2 = "#66a61e", ST8 = "#7570b3", ST13 = "#d95f02")

curve_dt <- rare_dt[!is.na(rarefied_observed)]
curve_dt <- merge(curve_dt, metrics[, .(sample, core, site_group, age_kyr, total_reads_from_counts, observed_taxa, qc_depth_class, low_reads_or_taxa, focus_age_band)], by = "sample", all.x = TRUE)

p_curve_all <- ggplot(curve_dt, aes(rarefaction_depth, rarefied_observed, group = sample)) +
  geom_line(aes(color = low_reads_or_taxa), alpha = 0.30, linewidth = 0.45) +
  scale_x_log10(labels = scales::label_number()) +
  scale_color_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey55"), na.value = "grey80") +
  labs(title = "Rarefaction curves by low-detection status", x = "Subsampled reads", y = "Expected observed taxa", color = "Low reads/taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "rarefaction_curves_low_detection.png"), p_curve_all, width = 8.5, height = 5.5, dpi = 160)

p_curve_core <- ggplot(curve_dt, aes(rarefaction_depth, rarefied_observed, group = sample, color = core)) +
  geom_line(alpha = 0.35, linewidth = 0.45) +
  facet_wrap(~ core) +
  scale_x_log10(labels = scales::label_number()) +
  scale_color_manual(values = core_cols) +
  labs(title = "Rarefaction curves by core", x = "Subsampled reads", y = "Expected observed taxa") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(OUT_FIG, "rarefaction_curves_by_core.png"), p_curve_core, width = 10, height = 6.5, dpi = 160)

p_focus <- ggplot(curve_dt[focus_age_band != "other"], aes(rarefaction_depth, rarefied_observed, group = sample, color = core)) +
  geom_line(alpha = 0.70, linewidth = 0.55) +
  facet_wrap(~ focus_age_band) +
  scale_x_log10(labels = scales::label_number()) +
  scale_color_manual(values = core_cols) +
  labs(title = "Rarefaction curves in focused age bands", x = "Subsampled reads", y = "Expected observed taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "rarefaction_curves_focus_age_bands.png"), p_focus, width = 9, height = 5, dpi = 160)

p_gap <- ggplot(metrics, aes(total_reads_from_counts, observed_taxa, color = qc_depth_class, shape = core)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_vline(xintercept = c(50000, 100000), linetype = "dashed", color = "grey45") +
  geom_hline(yintercept = c(150, 250), linetype = "dashed", color = "grey45") +
  scale_x_log10(labels = scales::label_number()) +
  labs(title = "Depth/richness QC classes", x = "Total reads", y = "Observed taxa", color = "QC class") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "depth_richness_qc_classes.png"), p_gap, width = 9, height = 5.8, dpi = 160)

p_slope <- ggplot(metrics, aes(age_kyr, rarefaction_slope_per_10k, color = core, shape = qc_depth_class)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_vline(xintercept = c(50, 100), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = core_cols) +
  labs(title = "Rarefaction slope across age", x = "Age (kyr)", y = "Added taxa per 10k additional reads near max depth") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "rarefaction_slope_by_age_core.png"), p_slope, width = 9, height = 5.5, dpi = 160)

p_class_core <- ggplot(core_class_summary, aes(core, pct_within_core, fill = qc_depth_class)) +
  geom_col(width = 0.75) +
  labs(title = "QC class composition by core", x = NULL, y = "Percent within core", fill = "QC class") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "qc_depth_class_by_core.png"), p_class_core, width = 8, height = 5, dpi = 160)

top_flags <- metrics[qc_depth_class %in% c("likely_undersequenced", "insufficient_for_call", "low_diversity_but_saturated"),
                     .(sample, core, age_kyr, total_reads_from_counts, observed_taxa, shannon_vegan, vegan_chao1,
                       chao1_gap, rarefaction_slope_per_10k, qc_depth_class)][order(qc_depth_class, core, age_kyr)]

sink(REPORT)
cat("# Rarefaction and Sequencing Depth QC Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Progress log: `InputQC/rarefaction_depth_qc/results/rarefaction_depth_qc.log`\n")
cat("- Packages: `phyloseq`, `vegan`, `iNEXT`\n\n")

cat("## Main interpretation\n\n")
cat("This analysis asks whether low-read / low-taxa samples appear undersequenced, or whether they look like genuinely low-diversity samples that have mostly saturated at their current depth.\n\n")
cat("Important caveat: Chao1 is weak for this processed count matrix because many samples have no singleton taxa; when Chao1 equals observed richness, the report relies more on rarefaction slope and observed depth/richness behavior.\n\n")

cat("## QC class summary\n\n")
cat("|QC class|n|pct|median reads|median observed taxa|median slope per 10k|\n")
cat("|---|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(class_summary))) {
  cat(sprintf("|%s|%d|%.1f|%.0f|%.0f|%.2f|\n",
              class_summary$qc_depth_class[i], class_summary$n[i], class_summary$pct[i],
              class_summary$median_reads[i], class_summary$median_observed_taxa[i],
              class_summary$median_slope_per_10k[i]))
}
cat("\n")

cat("## Threshold sensitivity summary\n\n")
cat("These thresholds should be treated as sensitivity candidates, not automatic production exclusions.\n\n")
cat("|threshold|flag value|n|pct|median reads|median observed taxa|pct likely undersequenced|pct low-diversity saturated|\n")
cat("|---|---|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(threshold_summary))) {
  cat(sprintf("|%s|%s|%d|%.1f|%.0f|%.0f|%.1f|%.1f|\n",
              threshold_summary$threshold_flag[i], as.character(threshold_summary$flag_value[i]),
              threshold_summary$n[i], threshold_summary$pct[i], threshold_summary$median_reads[i],
              threshold_summary$median_observed_taxa[i],
              threshold_summary$pct_likely_undersequenced[i],
              threshold_summary$pct_low_diversity_saturated[i]))
}
cat("\n")

cat("## Samples needing attention\n\n")
if (nrow(top_flags) == 0) {
  cat("No samples were assigned to the high-attention QC classes.\n\n")
} else {
  cat("|sample|core|age kyr|reads|observed taxa|slope per 10k|class|\n")
  cat("|---|---|---:|---:|---:|---:|---|\n")
  for (i in seq_len(min(80, nrow(top_flags)))) {
    cat(sprintf("|%s|%s|%.2f|%.0f|%.0f|%.2f|%s|\n",
                top_flags$sample[i], top_flags$core[i], top_flags$age_kyr[i],
                top_flags$total_reads_from_counts[i], top_flags$observed_taxa[i],
                top_flags$rarefaction_slope_per_10k[i], top_flags$qc_depth_class[i]))
  }
  if (nrow(top_flags) > 80) cat(sprintf("\nShowing first 80 of %d flagged samples. See the table output for all rows.\n\n", nrow(top_flags)))
}

cat("## CLR/log transform note\n\n")
cat("CLR and log transforms do not recover unobserved taxa. They can reduce scale effects among observed taxa, but if a sample lacks reads or detected taxa, a transform cannot infer the missing biological diversity. That is why rarefaction and depth QC must be evaluated before deciding whether a transformed matrix is trustworthy for network construction.\n\n")

cat("## Outputs\n\n")
cat("- Sample metrics: `InputQC/rarefaction_depth_qc/results/tables/sample_rarefaction_depth_metrics.tsv`\n")
cat("- Rarefaction checkpoints: `InputQC/rarefaction_depth_qc/results/tables/rarefaction_checkpoint_table.tsv`\n")
cat("- Threshold summary: `InputQC/rarefaction_depth_qc/results/tables/threshold_depth_summary.tsv`\n")
cat("- QC class summary: `InputQC/rarefaction_depth_qc/results/tables/qc_depth_class_summary.tsv`\n")
cat("- Figures: `InputQC/rarefaction_depth_qc/results/figures/`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
