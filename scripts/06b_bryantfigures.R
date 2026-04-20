#!/usr/bin/env Rscript
# 06b_bryantfigures.R — evaluation and manuscript figure

#load packages and config
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(here("config.R"))

# Path to pre-computed original analysis outputs used for new-vs-old comparisons.
# Update this if the original results are moved; it is intentionally not in
# config.R because it is system-specific and read-only reference data.
#OLD_WGCNA <- "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1"
FIGS      <- RESULTS$figures

theme_ms <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── Load new results ───────────────────────────────────────────────────────────

log_msg("Loading new results...")
meta      <- fread(UPSTREAM$metadata)
new_MEs   <- fread(file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv"))
new_mods  <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
new_pres  <- fread(file.path(RESULTS$stage1, "wgcna", "preservation.tsv"))
new_hmm   <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))
new_bic   <- fread(file.path(RESULTS$hmm, "bic_comparison.tsv"))
new_emp   <- fread(file.path(RESULTS$emp, "emp_sap_per_sample.tsv"))
new_tea   <- fread(file.path(RESULTS$tea, "tea_indices_per_sample.tsv"))
tea_corr  <- fread(file.path(RESULTS$tea, "tea_vs_emp_correlations.tsv"))
tea_clim  <- fread(file.path(RESULTS$tea, "tea_climate_models.tsv"))

# Merge metadata into new MEs (temp_complete loaded for potential future use)
new_trait <- merge(new_MEs, meta[, .(label, core, y_bp, mis, temp_complete)],
                   by.x = "sample", by.y = "label")
new_trait[, age_kyr := y_bp / 1000]
setnames(new_trait, "mis", "d18O")

# ── Load old results ───────────────────────────────────────────────────────────

# log_msg("Loading old results...")
# old_MEs  <- fread(file.path(OLD_WGCNA, "wgcna_prokaryotes", "module_eigengenes.tsv"))
# old_mods <- fread(file.path(OLD_WGCNA, "wgcna_prokaryotes", "module_assignments.tsv"))
# old_hmm  <- fread(file.path(OLD_WGCNA, "hmm_states", "hmm_states.tsv"))
# old_emp  <- fread(file.path(OLD_WGCNA, "emp_sap", "emp_sap_per_sample.tsv"))

# old_trait <- merge(old_MEs, meta[, .(label, core, y_bp, mis, temp_complete)],
#                    by.x = "sample", by.y = "label")
# old_trait[, age_kyr := y_bp / 1000]
# setnames(old_trait, "mis", "d18O")

# ── Fig 1: Module composition comparison ─────────────────────────────────────
# Melt module eigengene values for plotting
me_cols <- grep("^ME(?!grey)", names(new_trait), value = TRUE, perl = TRUE)

me_long <- melt(
  new_trait[, c("sample", "y_bp", "core", me_cols), with = FALSE],
  id.vars       = c("sample", "y_bp", "core"),
  variable.name = "module",
  value.name    = "eigengene"
)

me_long[, module := sub("^ME", "", module)]

p_me_vs_age <- ggplot(me_long, aes(x = y_bp, y = eigengene, colour = core)) +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~ module, scales = "free_y", nrow = 2) +
  scale_colour_brewer(palette = "Set2", name = "Core") +
  labs(
    title = "Module Eigen Taxa Loadings",
    x = "Years Before Present, BP",
    y = "Module Eigen Taxa Loadings"
  ) +
  theme_ms

