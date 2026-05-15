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
IN_DIR <- file.path(BASE, "results", "inputs")
OUT_TABLE <- file.path(BASE, "results", "tables")
OUT_FIG <- file.path(BASE, "results", "figures")
REPORT <- file.path(BASE, "INPUT_QC_DIAGNOSTIC_REPORT.md")
LOG <- file.path(BASE, "results", "input_qc_diagnostics.log")

dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) invisible(file.remove(LOG))

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

safe_cor <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

safe_lm_r2 <- function(y, x) {
  y <- as.numeric(y)
  ok <- is.finite(y) & complete.cases(x)
  if (sum(ok) < 5) return(NA_real_)
  x <- x[ok, , drop = FALSE]
  y <- y[ok]
  fit <- tryCatch(lm(y ~ ., data = as.data.frame(x)), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  summary(fit)$r.squared
}

variant_files <- c(
  current_taxon_centered_log = file.path(IN_DIR, "current_taxon_centered_log.rds"),
  deseq_length_log = file.path(IN_DIR, "deseq_length_log.rds"),
  sample_clr_raw = file.path(IN_DIR, "sample_clr_raw.rds"),
  log_depth_residualized = file.path(IN_DIR, "log_depth_residualized.rds")
)

missing <- variant_files[!file.exists(variant_files)]
if (length(missing)) stop("Missing input variant files: ", paste(names(missing), collapse = ", "))

log_msg("Loading metadata and depth summaries")
meta <- fread(file.path(IN_DIR, "sample_metadata_stage1.tsv"))
depth <- fread(file.path(IN_DIR, "sample_depth_summary.tsv"))
raw_counts <- readRDS(file.path(IN_DIR, "raw_counts.rds"))

meta_dt <- merge(meta, depth, by.x = "label", by.y = "sample", all.x = TRUE)
meta_dt[, `:=`(
  sample = label,
  age_kyr = y_bp / 1000,
  log_initial = log10(initial + 1),
  log_derep = log10(derep + 1),
  detected_taxa = as.integer(colSums(raw_counts > 0)[label]),
  shannon_raw = NA_real_
)]

count_by_sample <- t(raw_counts)
common_counts <- intersect(rownames(count_by_sample), meta_dt$sample)
shannon <- sapply(common_counts, function(s) {
  x <- count_by_sample[s, ]
  x <- x[x > 0]
  p <- x / sum(x)
  -sum(p * log(p))
})
meta_dt[match(common_counts, sample), shannon_raw := as.numeric(shannon)]

tech_covars <- c(
  "log_total_reads",
  "total_reads",
  "library_concentration",
  "log_initial",
  "log_derep",
  "avg_leng_initial",
  "avg_len_derep",
  "detected_taxa",
  "shannon_raw"
)
bio_covars <- c("age_kyr", "mis", "sst")
all_numeric_covars <- c(tech_covars, bio_covars)

sample_covariates <- meta_dt[, c(
  "sample", "core", "flowcell", "depth_in_core_cm",
  all_numeric_covars
), with = FALSE]
fwrite(sample_covariates, file.path(OUT_TABLE, "sample_technical_covariates.tsv"), sep = "\t")

log_msg("Computing covariate correlation table")
cov_cor <- rbindlist(lapply(all_numeric_covars, function(a) {
  rbindlist(lapply(all_numeric_covars, function(b) {
    data.table(
      covariate_a = a,
      covariate_b = b,
      pearson_r = safe_cor(meta_dt[[a]], meta_dt[[b]], "pearson"),
      spearman_rho = safe_cor(meta_dt[[a]], meta_dt[[b]], "spearman")
    )
  }))
}))
fwrite(cov_cor, file.path(OUT_TABLE, "technical_covariate_correlations.tsv"), sep = "\t")

log_msg("Computing PCA diagnostics for candidate inputs")
pc_assoc_rows <- list()
pc_score_rows <- list()
variance_rows <- list()
mean_sd_rows <- list()

for (variant in names(variant_files)) {
  log_msg("PCA diagnostics: ", variant)
  mat <- readRDS(variant_files[[variant]])
  common <- intersect(rownames(mat), meta_dt$sample)
  mat <- mat[common, , drop = FALSE]
  md <- meta_dt[match(common, sample)]

  keep <- apply(mat, 2, var, na.rm = TRUE) > 0
  mat <- mat[, keep, drop = FALSE]
  pc <- prcomp(mat, center = TRUE, scale. = FALSE)
  scores <- as.data.table(pc$x[, seq_len(min(5, ncol(pc$x))), drop = FALSE])
  scores[, sample := rownames(pc$x)]
  scores[, variant := variant]
  pc_score_rows[[length(pc_score_rows) + 1]] <- scores

  variance_rows[[length(variance_rows) + 1]] <- data.table(
    variant = variant,
    PC = paste0("PC", seq_len(min(10, length(pc$sdev)))),
    variance_explained = (pc$sdev^2 / sum(pc$sdev^2))[seq_len(min(10, length(pc$sdev)))]
  )

  pc_cols <- paste0("PC", seq_len(min(5, ncol(pc$x))))
  for (pc_col in pc_cols) {
    y <- scores[[pc_col]]
    for (covar in all_numeric_covars) {
      pc_assoc_rows[[length(pc_assoc_rows) + 1]] <- data.table(
        variant = variant,
        PC = pc_col,
        covariate = covar,
        covariate_class = fifelse(covar %in% tech_covars, "technical", "biological_proxy"),
        pearson_r = safe_cor(y, md[[covar]], "pearson"),
        spearman_rho = safe_cor(y, md[[covar]], "spearman")
      )
    }
    pc_assoc_rows[[length(pc_assoc_rows) + 1]] <- data.table(
      variant = variant,
      PC = pc_col,
      covariate = "core",
      covariate_class = "core_or_batch",
      pearson_r = sqrt(safe_lm_r2(y, model.matrix(~ core, md)[, -1, drop = FALSE])),
      spearman_rho = NA_real_
    )
  }

  feature_means <- colMeans(mat, na.rm = TRUE)
  feature_sds <- apply(mat, 2, sd, na.rm = TRUE)
  mean_sd_rows[[length(mean_sd_rows) + 1]] <- data.table(
    variant = variant,
    taxon = names(feature_means),
    feature_mean = feature_means,
    feature_sd = feature_sds,
    feature_cv_like = feature_sds / (abs(feature_means) + 1e-8)
  )
}

pc_assoc <- rbindlist(pc_assoc_rows, fill = TRUE)
pc_scores <- rbindlist(pc_score_rows, fill = TRUE)
variance_dt <- rbindlist(variance_rows, fill = TRUE)
mean_sd_dt <- rbindlist(mean_sd_rows, fill = TRUE)

pc_assoc[, abs_pearson := abs(pearson_r)]
fwrite(pc_assoc, file.path(OUT_TABLE, "input_pc_associations.tsv"), sep = "\t")
fwrite(pc_scores, file.path(OUT_TABLE, "input_pc_scores.tsv"), sep = "\t")
fwrite(variance_dt, file.path(OUT_TABLE, "input_pc_variance_explained.tsv"), sep = "\t")
fwrite(mean_sd_dt, file.path(OUT_TABLE, "input_feature_mean_sd.tsv"), sep = "\t")

top_assoc <- pc_assoc[PC %in% c("PC1", "PC2", "PC3")][
  order(variant, PC, -abs_pearson)
][, head(.SD, 5), by = .(variant, PC)]
fwrite(top_assoc, file.path(OUT_TABLE, "input_pc_top_associations.tsv"), sep = "\t")

log_msg("Writing figures")
heat_dt <- pc_assoc[PC %in% c("PC1", "PC2", "PC3") & covariate %in% all_numeric_covars]
heat_dt[, label := sprintf("%.2f", pearson_r)]
p_heat <- ggplot(heat_dt, aes(x = covariate, y = paste(variant, PC, sep = " / "), fill = pearson_r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(title = "Input PC associations with technical and biological covariates", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "input_pc_technical_heatmap.png"), p_heat, width = 12, height = 8, dpi = 160)

plot_dt <- merge(
  pc_scores[variant == "current_taxon_centered_log", .(sample, variant, PC1, PC2)],
  meta_dt[, .(sample, core, log_total_reads, age_kyr, mis, sst)],
  by = "sample",
  all.x = TRUE
)
p_ord_depth <- ggplot(plot_dt, aes(PC1, PC2, color = log_total_reads, shape = core)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_viridis_c(option = "magma") +
  labs(title = "Current input ordination colored by log total reads", color = "log10 reads") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "ordination_current_by_depth.png"), p_ord_depth, width = 7.5, height = 5.5, dpi = 160)

p_ord_core <- ggplot(plot_dt, aes(PC1, PC2, color = core)) +
  geom_point(size = 2.5, alpha = 0.9) +
  labs(title = "Current input ordination colored by core") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "ordination_current_by_core.png"), p_ord_core, width = 7.5, height = 5.5, dpi = 160)

p_var <- ggplot(variance_dt[PC %in% paste0("PC", 1:5)], aes(PC, variance_explained, group = variant, color = variant)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.7) +
  labs(title = "Variance explained by leading PCs", y = "Variance explained") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "input_pc_variance_explained.png"), p_var, width = 8, height = 5, dpi = 160)

assumption_register <- data.table(
  method_family = c(
    "Taxon-centered log current input",
    "Sample-wise CLR",
    "DESeq2 VST / rlog",
    "Microbiome size factors: GMPR/CSS/TMM",
    "Known covariate residualization",
    "SVA/RUV latent factor removal",
    "Filtering low-depth samples/taxa"
  ),
  useful_for = c(
    "Continuity with existing WGCNA and comparable module story",
    "Compositional sample centering and log-ratio geometry",
    "Mean-variance stabilization and size-factor normalization",
    "Zero-inflated microbiome-like count normalization",
    "Testing whether measured technical factors explain PCs",
    "Testing unknown unwanted variation",
    "Testing whether low-quality observations drive PC1"
  ),
  key_assumption_risk = c(
    "Does not remove sample-wide depth effects; can preserve technical axes",
    "Zeros and pseudocounts can dominate rare taxa; detection limits remain depth-dependent",
    "Built for RNA-seq-like negative-binomial count behavior; assumes size factors are meaningful for this system",
    "Often designed for microbiome differential abundance, not necessarily correlation-network inputs",
    "May erase real ancient-DNA preservation biology if depth/fragment metrics are biologically entangled",
    "Latent factors may represent true climate/core biology rather than nuisance variation",
    "Can improve quality but may bias age/core coverage and remove informative ancient samples"
  ),
  current_stance = c(
    "Anchor / current reference",
    "Required stress test, not automatically better",
    "Serious next candidate, but assumptions must be tracked",
    "Worth testing if packages are available; not first-line without diagnostics",
    "Control/sensitivity before production use",
    "Exploratory only until negative controls or clear nuisance factors are defined",
    "Likely important; must report sample/age/core loss"
  )
)
fwrite(assumption_register, file.path(OUT_TABLE, "method_assumption_register.tsv"), sep = "\t")

log_msg("Writing diagnostic report")
current_depth <- pc_assoc[variant == "current_taxon_centered_log" & PC == "PC1" & covariate == "log_total_reads", pearson_r]
worst_rows <- pc_assoc[PC %in% c("PC1", "PC2", "PC3")][order(-abs_pearson)][1:12]

sink(REPORT)
cat("# InputQC Diagnostic Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: diagnostic inventory before new normalization/correction choices\n")
cat("- Progress log: `InputQC/results/input_qc_diagnostics.log`\n\n")

cat("## Main finding\n\n")
cat(sprintf("The current input still has a strong PC1 association with log total reads: `r = %.3f`.\n\n", current_depth))
cat("This confirms that the next step should be diagnosis before correction. In ancient DNA, depth can be both technical and preservation-linked, so removing it blindly could erase real structure.\n\n")

cat("## Top PC associations\n\n")
cat("|variant|PC|covariate|class|Pearson r|Spearman rho|\n")
cat("|---|---|---|---|---:|---:|\n")
for (i in seq_len(nrow(worst_rows))) {
  r <- worst_rows[i]
  cat(sprintf("|%s|%s|%s|%s|%.3f|%.3f|\n",
              r$variant, r$PC, r$covariate, r$covariate_class,
              r$pearson_r, r$spearman_rho))
}
cat("\n")

cat("## Method assumption register\n\n")
cat("The goal is not to reject a method because it was designed elsewhere. The goal is to track what assumptions it may impose on ancient DNA data.\n\n")
cat("|method family|useful for|key assumption risk|current stance|\n")
cat("|---|---|---|---|\n")
for (i in seq_len(nrow(assumption_register))) {
  r <- assumption_register[i]
  cat(sprintf("|%s|%s|%s|%s|\n",
              r$method_family, r$useful_for, r$key_assumption_risk, r$current_stance))
}
cat("\n")

cat("## Recommended next move\n\n")
cat("Use these diagnostics to choose a small set of candidate corrections. The first serious candidates should be:\n\n")
cat("- DESeq2 VST with poscounts and reference-length normalization, as a count-model variance-stabilized input.\n")
cat("- sample-wise CLR with stricter filtering and pseudocount sensitivity, as the compositional baseline.\n")
cat("- depth/technical residualization as a control only, not as the default production input.\n")
cat("- sample/taxon filtering sensitivity, because low-depth samples and sparse taxa may drive PC1 as much as the transform does.\n\n")
cat("Do not promote any corrected input unless it reduces technical dominance while preserving climate/state/module biology.\n\n")

cat("## Outputs\n\n")
cat("- Covariates: `InputQC/results/tables/sample_technical_covariates.tsv`\n")
cat("- Covariate correlations: `InputQC/results/tables/technical_covariate_correlations.tsv`\n")
cat("- PC associations: `InputQC/results/tables/input_pc_associations.tsv`\n")
cat("- Top PC associations: `InputQC/results/tables/input_pc_top_associations.tsv`\n")
cat("- Method assumptions: `InputQC/results/tables/method_assumption_register.tsv`\n")
cat("- Heatmap: `InputQC/results/figures/input_pc_technical_heatmap.png`\n")
cat("- Ordination by depth: `InputQC/results/figures/ordination_current_by_depth.png`\n")
cat("- Ordination by core: `InputQC/results/figures/ordination_current_by_core.png`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
