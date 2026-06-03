#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

source(here("config.R"))
set.seed(PARAMS$seed)
options(stringsAsFactors = FALSE)

BASE <- here("InputQC")
IN_TABLE <- file.path(BASE, "results", "tables")
OUT_TABLE <- file.path(BASE, "results", "tables")
OUT_FIG <- file.path(BASE, "results", "figures")
REPORT <- file.path(BASE, "LOW_DETECTION_STRUCTURE_REPORT.md")
LOG <- file.path(BASE, "results", "low_detection_structure.log")

dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) invisible(file.remove(LOG))

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

safe_cor_test <- function(x, y, method = "spearman") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) {
    return(data.table(estimate = NA_real_, p_value = NA_real_, n = sum(ok)))
  }
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.table(estimate = unname(z$estimate), p_value = z$p.value, n = sum(ok))
}

safe_kw <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  if (sum(ok) < 4 || uniqueN(group[ok]) < 2) {
    return(data.table(statistic = NA_real_, p_value = NA_real_, n = sum(ok)))
  }
  z <- suppressWarnings(kruskal.test(value[ok] ~ as.factor(group[ok])))
  data.table(statistic = unname(z$statistic), p_value = z$p.value, n = sum(ok))
}

safe_wilcox <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  if (sum(ok) < 4 || uniqueN(group[ok]) != 2) {
    return(data.table(statistic = NA_real_, p_value = NA_real_, n = sum(ok)))
  }
  z <- suppressWarnings(wilcox.test(value[ok] ~ as.factor(group[ok]), exact = FALSE))
  data.table(statistic = unname(z$statistic), p_value = z$p.value, n = sum(ok))
}

num_min <- function(x) as.numeric(min(x, na.rm = TRUE))
num_median <- function(x) as.numeric(median(x, na.rm = TRUE))
num_max <- function(x) as.numeric(max(x, na.rm = TRUE))

log_msg("Loading sample covariates and PC scores")
cov <- fread(file.path(IN_TABLE, "sample_technical_covariates.tsv"))
pc <- fread(file.path(IN_TABLE, "input_pc_scores.tsv"))[
  variant == "current_taxon_centered_log",
  .(sample, PC1, PC2, PC3)
]
dt <- merge(cov, pc, by = "sample", all.x = TRUE)

dt[, site_group := fifelse(core %in% c("GeoB25202_R1", "GeoB25202_R2"), "GeoB25202", core)]

# d18O is stored as `mis`; use a simple project convention for descriptive glacial classes.
dt[, glacial_class := fifelse(is.na(mis), NA_character_,
                              fifelse(mis >= 4.0, "glacial_like", "interglacial_like"))]
dt[, dedup_fraction := fifelse(is.finite(log_initial) & is.finite(log_derep),
                               (10^log_derep - 1) / (10^log_initial - 1),
                               NA_real_)]

read_q10 <- quantile(dt$total_reads, 0.10, na.rm = TRUE, names = FALSE)
read_q20 <- quantile(dt$total_reads, 0.20, na.rm = TRUE, names = FALSE)
taxa_q10 <- quantile(dt$detected_taxa, 0.10, na.rm = TRUE, names = FALSE)
taxa_q20 <- quantile(dt$detected_taxa, 0.20, na.rm = TRUE, names = FALSE)

dt[, `:=`(
  low_reads_50k = total_reads < 50000,
  low_reads_100k = total_reads < 100000,
  low_reads_q10 = total_reads <= read_q10,
  low_reads_q20 = total_reads <= read_q20,
  low_taxa_150 = detected_taxa < 150,
  low_taxa_250 = detected_taxa < 250,
  low_taxa_q10 = detected_taxa <= taxa_q10,
  low_taxa_q20 = detected_taxa <= taxa_q20,
  low_reads_or_taxa = total_reads < 100000 | detected_taxa < 250,
  low_reads_and_taxa = total_reads < 100000 & detected_taxa < 250
)]

fwrite(dt, file.path(OUT_TABLE, "low_detection_sample_table.tsv"), sep = "\t")

group_vars <- c(
  "low_reads_50k", "low_reads_100k", "low_reads_q10", "low_reads_q20",
  "low_taxa_150", "low_taxa_250", "low_taxa_q10", "low_taxa_q20",
  "low_reads_or_taxa", "low_reads_and_taxa"
)

summarize_group <- function(group_var) {
  gv <- group_var
  out <- dt[, .(
    n = .N,
    pct = .N / nrow(dt) * 100,
    age_min = num_min(age_kyr),
    age_median = num_median(age_kyr),
    age_max = num_max(age_kyr),
    d18O_median = num_median(mis),
    sst_median = num_median(sst),
    total_reads_median = num_median(total_reads),
    detected_taxa_median = num_median(detected_taxa),
    shannon_median = num_median(shannon_raw),
    library_concentration_median = num_median(library_concentration),
    avg_len_derep_median = num_median(avg_len_derep),
    dedup_fraction_median = num_median(dedup_fraction),
    PC1_median = num_median(PC1),
    PC2_median = num_median(PC2)
  ), by = .(group = get(gv))][order(group)]
  out[, group_var := gv]
  setcolorder(out, c("group_var", "group"))
  out
}

group_summary <- rbindlist(lapply(group_vars, summarize_group), fill = TRUE)
fwrite(group_summary, file.path(OUT_TABLE, "low_detection_group_summary.tsv"), sep = "\t")

core_summary <- rbindlist(lapply(group_vars, function(gv) {
  out <- dt[, .(n = .N), by = .(group = get(gv), core)]
  out[, group_var := gv]
  out[, pct_within_group := n / sum(n) * 100, by = .(group_var, group)]
  setcolorder(out, c("group_var", "group", "core"))
  out
}), fill = TRUE)
fwrite(core_summary, file.path(OUT_TABLE, "low_detection_core_enrichment.tsv"), sep = "\t")

site_summary <- rbindlist(lapply(group_vars, function(gv) {
  out <- dt[, .(n = .N), by = .(group = get(gv), site_group, glacial_class)]
  out[, group_var := gv]
  out[, pct_within_group := n / sum(n) * 100, by = .(group_var, group)]
  setcolorder(out, c("group_var", "group", "site_group", "glacial_class"))
  out
}), fill = TRUE)
fwrite(site_summary, file.path(OUT_TABLE, "low_detection_site_glacial_summary.tsv"), sep = "\t")

metrics <- c(
  "detected_taxa", "total_reads", "log_total_reads", "shannon_raw",
  "library_concentration", "avg_leng_initial", "avg_len_derep",
  "dedup_fraction", "PC1", "PC2"
)
continuous_predictors <- c("age_kyr", "mis", "sst", "depth_in_core_cm")

assoc_cor <- rbindlist(lapply(metrics, function(metric) {
  rbindlist(lapply(continuous_predictors, function(pred) {
    out <- safe_cor_test(dt[[metric]], dt[[pred]], "spearman")
    out[, `:=`(test = "spearman", response = metric, predictor = pred)]
    out
  }))
}), fill = TRUE)

assoc_group <- rbindlist(lapply(metrics, function(metric) {
  kw_core <- safe_kw(dt[[metric]], dt$core)
  kw_core[, `:=`(test = "kruskal_core", response = metric, predictor = "core")]
  kw_site <- safe_kw(dt[[metric]], dt$site_group)
  kw_site[, `:=`(test = "kruskal_site_group", response = metric, predictor = "site_group")]
  wx_glacial <- safe_wilcox(dt[[metric]], dt$glacial_class)
  wx_glacial[, `:=`(test = "wilcox_glacial_class", response = metric, predictor = "glacial_class")]
  rbindlist(list(kw_core, kw_site, wx_glacial), fill = TRUE)
}), fill = TRUE)

assoc_tests <- rbindlist(list(assoc_cor, assoc_group), fill = TRUE)
assoc_tests[, fdr := p.adjust(p_value, method = "BH")]
setcolorder(assoc_tests, c("test", "response", "predictor", "estimate", "statistic", "p_value", "fdr", "n"))
fwrite(assoc_tests, file.path(OUT_TABLE, "low_detection_association_tests.tsv"), sep = "\t")

log_msg("Writing figures")
core_cols <- c(
  GeoB25202_R1 = "#1b9e77",
  GeoB25202_R2 = "#66a61e",
  ST8 = "#7570b3",
  ST13 = "#d95f02"
)

p_detect_age <- ggplot(dt, aes(age_kyr, detected_taxa, color = core)) +
  geom_point(aes(shape = glacial_class), size = 2.3, alpha = 0.85) +
  geom_hline(yintercept = c(150, 250), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = core_cols) +
  labs(title = "Detected taxa across age", x = "Age (kyr)", y = "Detected taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "detected_taxa_by_age_core.png"), p_detect_age, width = 8, height = 5, dpi = 160)

p_reads_age <- ggplot(dt, aes(age_kyr, log_total_reads, color = core)) +
  geom_point(aes(shape = glacial_class), size = 2.3, alpha = 0.85) +
  geom_hline(yintercept = log10(c(50000, 100000)), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = core_cols) +
  labs(title = "Read depth across age", x = "Age (kyr)", y = "log10 total reads") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "reads_by_age_core.png"), p_reads_age, width = 8, height = 5, dpi = 160)

p_shannon_age <- ggplot(dt, aes(age_kyr, shannon_raw, color = core)) +
  geom_point(aes(shape = glacial_class), size = 2.3, alpha = 0.85) +
  scale_color_manual(values = core_cols) +
  labs(title = "Raw Shannon diversity across age", x = "Age (kyr)", y = "Shannon diversity") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "shannon_by_age_core.png"), p_shannon_age, width = 8, height = 5, dpi = 160)

p_detect_d18o <- ggplot(dt, aes(mis, detected_taxa, color = core)) +
  geom_point(aes(shape = glacial_class), size = 2.3, alpha = 0.85) +
  geom_vline(xintercept = 4.0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(150, 250), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = core_cols) +
  labs(title = "Detected taxa vs d18O", x = "d18O proxy", y = "Detected taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "detected_taxa_by_d18O_core.png"), p_detect_d18o, width = 8, height = 5, dpi = 160)

p_detect_reads <- ggplot(dt, aes(log_total_reads, detected_taxa, color = site_group)) +
  geom_point(aes(shape = glacial_class), size = 2.3, alpha = 0.85) +
  geom_vline(xintercept = log10(c(50000, 100000)), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(150, 250), linetype = "dashed", color = "grey40") +
  labs(title = "Detected taxa vs read depth", x = "log10 total reads", y = "Detected taxa") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "detected_vs_reads_by_core_glacial.png"), p_detect_reads, width = 8, height = 5, dpi = 160)

pc_long <- melt(dt[, .(sample, core, site_group, glacial_class, PC1, detected_taxa, log_total_reads, shannon_raw)],
                id.vars = c("sample", "core", "site_group", "glacial_class", "PC1"),
                variable.name = "metric", value.name = "value")
p_pc1 <- ggplot(pc_long, aes(value, PC1, color = core)) +
  geom_point(aes(shape = glacial_class), size = 2.2, alpha = 0.85) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_color_manual(values = core_cols) +
  labs(title = "PC1 relationship to detection metrics", x = NULL, y = "Current input PC1") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "pc1_vs_detection_metrics.png"), p_pc1, width = 10, height = 4.8, dpi = 160)

box_dt <- melt(dt[, .(sample, core, site_group, glacial_class, detected_taxa, log_total_reads, shannon_raw, dedup_fraction, avg_len_derep)],
               id.vars = c("sample", "core", "site_group", "glacial_class"),
               variable.name = "metric", value.name = "value")
p_box <- ggplot(box_dt, aes(core, value, fill = glacial_class)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.75) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c(glacial_like = "#6baed6", interglacial_like = "#fd8d3c"), na.value = "grey70") +
  labs(title = "QC metrics by core and glacial class", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_FIG, "qc_metrics_by_core_glacial.png"), p_box, width = 10, height = 7, dpi = 160)

log_msg("Classifying low-detection structure")
primary <- dt[low_reads_or_taxa == TRUE]
primary_core <- primary[, .N, by = core][order(-N)]
primary_site <- primary[, .N, by = site_group][order(-N)]
primary_glacial <- primary[, .N, by = glacial_class][order(-N)]
primary_age_min <- primary[, min(age_kyr, na.rm = TRUE)]
primary_age_max <- primary[, max(age_kyr, na.rm = TRUE)]
primary_pct_st8 <- primary[core == "ST8", .N] / nrow(primary)
primary_pct_glacial <- primary[glacial_class == "glacial_like", .N] / nrow(primary)

classification <- "mixed"
if (primary_pct_st8 > 0.75) classification <- "mostly core/site structured"
if (primary_pct_glacial > 0.75) classification <- "mostly glacial/paleoenvironment structured"
if (primary_pct_st8 <= 0.75 && primary_pct_glacial <= 0.75 &&
    nrow(primary_core) >= 3 && primary_age_min < 20 && primary_age_max > 100) {
  classification <- "mixed technical and paleoenvironmental/core structure"
}

sink(REPORT)
cat("# Low-Detection Structure Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: qualitative/descriptive investigation of the convergent low-depth/low-detection ordination region\n")
cat("- Progress log: `InputQC/results/low_detection_structure.log`\n\n")

cat("## Main interpretation\n\n")
cat("The convergent region in the current ordination is best treated as a low-detection axis, not as a PC-threshold removal target.\n\n")
cat(sprintf("Using the primary descriptive group `total_reads < 100000 OR detected_taxa < 250`, the low-detection group contains `%d` samples across ages %.2f-%.2f kyr.\n\n",
            nrow(primary), primary_age_min, primary_age_max))
cat(sprintf("Current qualitative classification: **%s**.\n\n", classification))
cat("This means filtering may be useful as a sensitivity analysis, but it is not yet justified as a production exclusion rule. The samples may include technical weakness, preservation structure, and real paleoenvironmental signal.\n\n")

cat("## Primary group composition\n\n")
cat("Core composition for `total_reads < 100000 OR detected_taxa < 250`:\n\n")
cat("|core|n|\n|---|---:|\n")
for (i in seq_len(nrow(primary_core))) {
  cat(sprintf("|%s|%d|\n", primary_core$core[i], primary_core$N[i]))
}
cat("\nSite composition:\n\n")
cat("|site_group|n|\n|---|---:|\n")
for (i in seq_len(nrow(primary_site))) {
  cat(sprintf("|%s|%d|\n", primary_site$site_group[i], primary_site$N[i]))
}
cat("\nGlacial class composition:\n\n")
cat("|glacial_class|n|\n|---|---:|\n")
for (i in seq_len(nrow(primary_glacial))) {
  cat(sprintf("|%s|%d|\n", primary_glacial$glacial_class[i], primary_glacial$N[i]))
}
cat("\n")

cat("## Filter sensitivity framing\n\n")
cat("Reasonable candidate filters to test later, not apply immediately:\n\n")
cat("- `detected_taxa >= 150`: removes the most extreme low-richness tail.\n")
cat("- `detected_taxa >= 250`: stronger richness filter, close to the visual PC1>50 group but based on independent QC.\n")
cat("- `total_reads >= 50000`: conservative depth filter.\n")
cat("- `total_reads >= 100000`: stronger depth filter.\n")
cat("- combined filter `total_reads >= 50000 AND detected_taxa >= 150`: balanced first sensitivity candidate.\n\n")
cat("A PC1 threshold should remain diagnostic only. It should not be used as the filtering rule.\n\n")

cat("## Association tests\n\n")
cat("Association tests are descriptive because samples are time-ordered and cores are not independent. Use them to guide plots and sensitivity analyses, not as final inference.\n\n")
cat("See `InputQC/results/tables/low_detection_association_tests.tsv`.\n\n")

cat("## Outputs\n\n")
cat("- Sample table: `InputQC/results/tables/low_detection_sample_table.tsv`\n")
cat("- Group summaries: `InputQC/results/tables/low_detection_group_summary.tsv`\n")
cat("- Core enrichment: `InputQC/results/tables/low_detection_core_enrichment.tsv`\n")
cat("- Site/glacial summary: `InputQC/results/tables/low_detection_site_glacial_summary.tsv`\n")
cat("- Association tests: `InputQC/results/tables/low_detection_association_tests.tsv`\n")
cat("- Figures: `InputQC/results/figures/*low*`, `*detected*`, `*reads*`, `*shannon*`, `*pc1*`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
