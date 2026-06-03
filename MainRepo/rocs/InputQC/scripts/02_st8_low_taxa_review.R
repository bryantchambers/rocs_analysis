#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

source(here("config.R"))
set.seed(PARAMS$seed)
options(stringsAsFactors = FALSE)

BASE <- here("InputQC", "st8_low_taxa_review")
OUT_TABLE <- file.path(BASE, "results", "tables")
OUT_FIG <- file.path(BASE, "results", "figures")
REPORT <- file.path(BASE, "ST8_LOW_TAXA_REVIEW.md")
LOG <- file.path(BASE, "results", "st8_low_taxa_review.log")

dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) invisible(file.remove(LOG))

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

summary_row <- function(x, label) {
  data.table(
    group = label,
    n = nrow(x),
    median_library_concentration = median(x$library_concentration, na.rm = TRUE),
    p25_library_concentration = quantile(x$library_concentration, 0.25, na.rm = TRUE, names = FALSE),
    p75_library_concentration = quantile(x$library_concentration, 0.75, na.rm = TRUE, names = FALSE),
    median_total_reads = median(x$total_reads, na.rm = TRUE),
    median_detected_taxa = median(x$detected_taxa, na.rm = TRUE),
    median_shannon = median(x$shannon_raw, na.rm = TRUE),
    median_age_kyr = median(x$age_kyr, na.rm = TRUE)
  )
}

log_msg("Loading InputQC tables and raw counts")
raw_counts <- readRDS(here("InputQC", "results", "inputs", "raw_counts.rds"))
low <- fread(here("InputQC", "results", "tables", "low_detection_sample_table.tsv"))
rare <- fread(here("InputQC", "rarefaction_depth_qc", "results", "tables", "sample_rarefaction_depth_metrics.tsv"))
tax_meta <- fread(here("InputQC", "results", "inputs", "prokaryotes_taxa_metadata.tsv"))

dt <- merge(
  low,
  rare[, .(sample, observed_taxa, shannon_vegan, qc_depth_class, rarefaction_slope_per_10k)],
  by = "sample",
  all.x = TRUE
)

st8 <- dt[core == "ST8"]
st8_low <- st8[low_taxa_250 %in% TRUE]
st8_not_low <- st8[!(low_taxa_250 %in% TRUE)]
non_st8 <- dt[core != "ST8"]

group_summary <- rbindlist(list(
  summary_row(st8_low, "ST8 low taxa (<250)"),
  summary_row(st8_not_low, "ST8 not low taxa"),
  summary_row(non_st8, "non-ST8"),
  summary_row(dt, "all samples")
), fill = TRUE)

sample_table <- st8_low[order(age_kyr), .(
  sample, age_kyr, depth_in_core_cm, total_reads, detected_taxa,
  observed_taxa, shannon_raw, library_concentration, log_initial,
  log_derep, avg_leng_initial, avg_len_derep, qc_depth_class,
  rarefaction_slope_per_10k
)]

samples <- intersect(sample_table$sample, colnames(raw_counts))
sub <- raw_counts[, samples, drop = FALSE]

taxa_summary <- data.table(
  taxon = rownames(sub),
  samples_present = rowSums(sub > 0),
  total_counts = rowSums(sub),
  max_count = apply(sub, 1, max),
  mean_count_when_present = apply(sub, 1, function(x) mean(x[x > 0]))
)[samples_present > 0][order(-samples_present, -total_counts)]
taxa_summary[, pct_low_samples_present := samples_present / length(samples) * 100]
taxa_summary <- merge(taxa_summary, tax_meta, by = "taxon", all.x = TRUE)
setorder(taxa_summary, -samples_present, -total_counts)

per_sample_taxa <- rbindlist(lapply(samples, function(s) {
  x <- sub[, s]
  out <- data.table(sample = s, taxon = rownames(sub), count = as.numeric(x))[count > 0]
  out[, relative_abundance := count / sum(count)]
  out <- merge(out, tax_meta, by = "taxon", all.x = TRUE)
  out[order(-count)]
}), fill = TRUE)

top_per_sample <- per_sample_taxa[, head(.SD, 15), by = sample]
setorder(top_per_sample, sample, -count)

read_distribution <- rbindlist(lapply(samples, function(s) {
  x <- sort(as.numeric(sub[, s][sub[, s] > 0]), decreasing = TRUE)
  total <- sum(x)
  p <- x / total
  cum <- cumsum(p)
  data.table(
    sample = s,
    detected_taxa = length(x),
    total_reads = total,
    top1_pct = 100 * p[1],
    top5_pct = 100 * sum(head(p, 5)),
    top10_pct = 100 * sum(head(p, 10)),
    taxa_for_50pct = which(cum >= 0.50)[1],
    taxa_for_80pct = which(cum >= 0.80)[1],
    taxa_for_90pct = which(cum >= 0.90)[1],
    pielou_evenness = (-sum(p * log(p))) / log(length(p)),
    gini = if (length(x) > 1) {
      n <- length(x)
      sx <- sort(x)
      (2 * sum(seq_len(n) * sx) / (n * sum(sx))) - (n + 1) / n
    } else {
      0
    }
  )
}), fill = TRUE)
read_distribution <- merge(
  sample_table[, .(sample, age_kyr, library_concentration, qc_depth_class)],
  read_distribution,
  by = "sample",
  all.x = TRUE
)
setorder(read_distribution, age_kyr)

taxa_function_summary <- taxa_summary[, .(
  taxa_n = .N,
  median_samples_present = median(samples_present, na.rm = TRUE),
  total_counts = sum(total_counts, na.rm = TRUE)
), by = .(domain, phylum, functional_group, signal_source)][order(-total_counts)]

fwrite(group_summary, file.path(OUT_TABLE, "st8_low_taxa_group_summary.tsv"), sep = "\t")
fwrite(sample_table, file.path(OUT_TABLE, "st8_low_taxa_sample_table.tsv"), sep = "\t")
fwrite(taxa_summary, file.path(OUT_TABLE, "st8_low_taxa_taxa_summary.tsv"), sep = "\t")
fwrite(top_per_sample, file.path(OUT_TABLE, "st8_low_taxa_top_taxa_per_sample.tsv"), sep = "\t")
fwrite(taxa_function_summary, file.path(OUT_TABLE, "st8_low_taxa_function_summary.tsv"), sep = "\t")
fwrite(read_distribution, file.path(OUT_TABLE, "st8_low_taxa_read_distribution.tsv"), sep = "\t")

log_msg("Writing figures")
plot_dt <- copy(dt)
plot_dt[, st8_low_taxa_group := fifelse(core == "ST8" & low_taxa_250 %in% TRUE, "ST8 low taxa",
                                 fifelse(core == "ST8", "ST8 not low taxa", "non-ST8"))]

p_conc <- ggplot(plot_dt, aes(st8_low_taxa_group, library_concentration, fill = st8_low_taxa_group)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.5) +
  labs(title = "Library concentration in ST8 low-taxa samples", x = NULL, y = "Library concentration") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(OUT_FIG, "library_concentration_st8_low_taxa.png"), p_conc, width = 7, height = 5, dpi = 160)

p_age <- ggplot(st8, aes(age_kyr, library_concentration, color = low_taxa_250)) +
  geom_point(aes(size = detected_taxa), alpha = 0.85) +
  geom_vline(xintercept = c(50, 100), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey45"), na.value = "grey80") +
  labs(title = "ST8 library concentration across age", x = "Age (kyr)", y = "Library concentration", color = "Low taxa <250", size = "Detected taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "st8_library_concentration_by_age.png"), p_age, width = 8, height = 5, dpi = 160)

top_taxa <- taxa_summary[1:min(25, .N)]
p_taxa <- ggplot(top_taxa, aes(reorder(taxon, samples_present), samples_present)) +
  geom_col(fill = "#4c78a8") +
  coord_flip() +
  labs(title = "Most recurrent taxa in ST8 low-taxa samples", x = NULL, y = "Samples present") +
  theme_minimal(base_size = 9)
ggsave(file.path(OUT_FIG, "recurrent_taxa_st8_low_taxa.png"), p_taxa, width = 8, height = 7, dpi = 160)

top_heat_taxa <- taxa_summary[1:min(40, .N), taxon]
heat_dt <- per_sample_taxa[taxon %in% top_heat_taxa]
heat_dt <- merge(heat_dt, sample_table[, .(sample, age_kyr)], by = "sample", all.x = TRUE)
heat_dt[, sample_label := sprintf("%s\n%.1f kya", sample, age_kyr)]
p_heat <- ggplot(heat_dt, aes(sample_label, taxon, fill = log10(count + 1))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Top recurrent taxa across ST8 low-taxa samples", x = NULL, y = NULL, fill = "log10(count+1)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "st8_low_taxa_top_taxa_heatmap.png"), p_heat, width = 12, height = 8, dpi = 160)

dist_long <- melt(
  read_distribution,
  id.vars = c("sample", "age_kyr", "detected_taxa", "library_concentration"),
  measure.vars = c("top1_pct", "top5_pct", "top10_pct"),
  variable.name = "read_share_metric",
  value.name = "pct_reads"
)
p_dist <- ggplot(dist_long, aes(age_kyr, pct_reads, color = read_share_metric)) +
  geom_point(aes(size = detected_taxa), alpha = 0.85) +
  geom_line(alpha = 0.55) +
  geom_vline(xintercept = c(50, 100), linetype = "dashed", color = "grey55") +
  labs(title = "Read concentration across detected taxa in ST8 low-taxa samples",
       x = "Age (kyr)", y = "Percent of reads", color = "Metric", size = "Detected taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "st8_low_taxa_read_concentration_by_age.png"), p_dist, width = 8.5, height = 5.2, dpi = 160)

p_even <- ggplot(read_distribution, aes(age_kyr, pielou_evenness, color = library_concentration)) +
  geom_point(aes(size = detected_taxa), alpha = 0.9) +
  geom_vline(xintercept = c(50, 100), linetype = "dashed", color = "grey55") +
  scale_color_viridis_c(option = "magma") +
  labs(title = "Evenness in ST8 low-taxa samples",
       x = "Age (kyr)", y = "Pielou evenness", color = "Library concentration", size = "Detected taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "st8_low_taxa_evenness_by_age.png"), p_even, width = 8.5, height = 5.2, dpi = 160)

sink(REPORT)
cat("# ST8 Low-Taxa Sample Review\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Definition: ST8 samples with `detected_taxa < 250`\n")
cat("- Progress log: `InputQC/st8_low_taxa_review/results/st8_low_taxa_review.log`\n\n")

cat("## Main findings\n\n")
cat(sprintf("There are `%d` ST8 low-taxa samples under the `<250 detected taxa` threshold.\n\n", nrow(sample_table)))
cat("Their library concentration is substantially lower than the rest of ST8 and the overall dataset, but not uniformly low; several low-taxa samples still have moderate library concentration. That argues for a mixed process: DNA/library yield contributes strongly, but it probably does not explain every low-taxa ST8 sample by itself.\n\n")

cat("## DNA/library concentration comparison\n\n")
cat("|group|n|median library concentration|IQR|median reads|median detected taxa|median Shannon|median age kyr|\n")
cat("|---|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(group_summary))) {
  cat(sprintf("|%s|%d|%.2f|%.2f-%.2f|%.0f|%.0f|%.2f|%.2f|\n",
              group_summary$group[i], group_summary$n[i],
              group_summary$median_library_concentration[i],
              group_summary$p25_library_concentration[i],
              group_summary$p75_library_concentration[i],
              group_summary$median_total_reads[i],
              group_summary$median_detected_taxa[i],
              group_summary$median_shannon[i],
              group_summary$median_age_kyr[i]))
}
cat("\n")

cat("## Most recurrent taxa\n\n")
cat("|taxon|samples present|pct samples|total counts|domain|phylum|functional group|signal source|\n")
cat("|---|---:|---:|---:|---|---|---|---|\n")
for (i in seq_len(min(25, nrow(taxa_summary)))) {
  cat(sprintf("|%s|%d|%.1f|%.0f|%s|%s|%s|%s|\n",
              taxa_summary$taxon[i], taxa_summary$samples_present[i],
              taxa_summary$pct_low_samples_present[i], taxa_summary$total_counts[i],
              taxa_summary$domain[i], taxa_summary$phylum[i],
              taxa_summary$functional_group[i], taxa_summary$signal_source[i]))
}
cat("\n")

dist_summary <- read_distribution[, .(
  median_top1_pct = median(top1_pct, na.rm = TRUE),
  median_top5_pct = median(top5_pct, na.rm = TRUE),
  median_top10_pct = median(top10_pct, na.rm = TRUE),
  median_taxa_for_50pct = median(taxa_for_50pct, na.rm = TRUE),
  median_taxa_for_80pct = median(taxa_for_80pct, na.rm = TRUE),
  median_taxa_for_90pct = median(taxa_for_90pct, na.rm = TRUE),
  median_pielou_evenness = median(pielou_evenness, na.rm = TRUE),
  median_gini = median(gini, na.rm = TRUE)
)]

cat("## Read distribution across detected taxa\n\n")
cat(sprintf("Across ST8 low-taxa samples, the median top taxon contains `%.1f%%` of reads, the median top 5 taxa contain `%.1f%%`, and the median top 10 taxa contain `%.1f%%`.\n\n",
            dist_summary$median_top1_pct, dist_summary$median_top5_pct, dist_summary$median_top10_pct))
cat(sprintf("The median sample needs `%d` taxa to explain 50%% of reads, `%d` taxa to explain 80%%, and `%d` taxa to explain 90%%. Median Pielou evenness is `%.2f`, so most samples are not collapsed into a single dominant taxon, although several sparse samples are strongly concentrated in their top few taxa.\n\n",
            as.integer(dist_summary$median_taxa_for_50pct),
            as.integer(dist_summary$median_taxa_for_80pct),
            as.integer(dist_summary$median_taxa_for_90pct),
            dist_summary$median_pielou_evenness))

cat("## Outputs\n\n")
cat("- Sample table: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_sample_table.tsv`\n")
cat("- Taxa summary: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_taxa_summary.tsv`\n")
cat("- Top taxa per sample: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_top_taxa_per_sample.tsv`\n")
cat("- Function summary: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_function_summary.tsv`\n")
cat("- Read distribution: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_read_distribution.tsv`\n")
cat("- Figures: `InputQC/st8_low_taxa_review/results/figures/`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
