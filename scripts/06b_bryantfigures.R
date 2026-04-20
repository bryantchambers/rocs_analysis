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

# ── Module relative abundance by sample ──────────────────────────────────────

log_msg("Computing module relative abundances...")

# ── 1. Load raw reads and filter ──────────────────────────────────────────────

prok_raw <- fread(UPSTREAM$tax_damage)

prok <- prok_raw[
  is_dmg == "Damaged" &
  domain %in% c("d__Archaea", "d__Bacteria", "Viruses") &
  label %in% meta$label
]

# Aggregate by subspecies × sample
prok_agg <- prok[, .(
  n_reads = sum(n_reads)
), by = .(subspecies, label)]

log_msg(sprintf("  Raw: %d taxa × %d samples", 
                uniqueN(prok_agg$subspecies), uniqueN(prok_agg$label)))

# ── 2. Load and merge metadata ───────────────────────────────────────────────

# Module assignments
mods <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
setnames(mods, "taxon", "subspecies")

# Merge modules into aggregated reads
prok_mod <- merge(prok_agg, mods, by = "subspecies", all.x = TRUE)
prok_mod[is.na(module), module := "grey"]

log_msg(sprintf("  After module assignment: %d records", nrow(prok_mod)))

# ── 3. Normalize to relative abundance per sample ────────────────────────────

# Total reads per sample
prok_mod[, total_reads := sum(n_reads), by = label]
prok_mod[, rel_abund := n_reads / total_reads]

log_msg("  Relative abundance computed (sum per sample = 1.0)")

# ── 4. Aggregate by module × sample ──────────────────────────────────────────

module_abund <- prok_mod[, .(
  n_reads = sum(n_reads),
  rel_abund = sum(rel_abund)
), by = .(module, label)]

# Merge metadata (y_bp, core)
module_abund <- merge(
  module_abund,
  meta[, .(label, core, y_bp, age_kyr = y_bp / 1000)],
  by = "label"
)

# Verify relative abundance per sample sums to 1 (within floating-point tolerance)
abund_check <- module_abund[, .(total_rel = sum(rel_abund)), by = label]
log_msg(sprintf("  Abundance check: min=%.6f, max=%.6f (should be ~1.0)",
                min(abund_check$total_rel), max(abund_check$total_rel)))

# ── 5. Order samples by age and cores ────────────────────────────────────────

module_abund <- module_abund[order(core, -y_bp)]
module_abund[, sample_order := paste0(core, "_", sprintf("%05d", rank(-y_bp)), by = core)]

# ── 6. Stacked bar plot: Module composition across cores ────────────────────

# Define module colors (consistent with theme)
module_colors_pal <- c(
  turquoise = "#1B9E77",
  blue      = "#377EB8",
  brown     = "#A65628",
  yellow    = "#E6AB02",
  green     = "#66A61E",
  red       = "#E7298A",
  grey      = "#CCCCCC"
)

# Optional: bin by age for smoother visualization
module_abund[, age_bin := floor(age_kyr / 10) * 10 + 5]
module_abund_binned <- module_abund[, .(
  rel_abund = sum(rel_abund)
), by = .(module, core, age_bin)]
# Then plot using age_bin on x-axis

# p_module_abund <- ggplot(module_abund[module != "grey"], 
#                           aes(x = reorder(label, -y_bp), y = rel_abund, fill = module)) +
#   geom_col(colour = NA, width = 0.85) +
#   facet_wrap(~ core, scales = "free_x", nrow = 1) +
#   scale_fill_manual(values = module_colors_pal, name = "Module") +
#   labs(
#     title = "Module relative abundance by sample",
#     x = "Sample (ordered by age)",
#     y = "Relative abundance"
#   ) +
#   theme_ms +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#     strip.text = element_text(size = 9, face = "bold")
#   )

log_msg("Module abundance plot created")

# ── Area plot of module relative abundance vs y_bp ───────────────────────────

module_abund <- module_abund[order(core, y_bp)]

module_colors_pal <- c(
  turquoise = "#1B9E77",
  blue      = "#377EB8",
  brown     = "#A65628",
  yellow    = "#E6AB02",
  green     = "#66A61E",
  red       = "#E7298A",
  grey      = "#CCCCCC"
)

p_module_abund <- ggplot(
  module_abund[module != "grey"],
  aes(x = y_bp, y = rel_abund, fill = module, group = module)
) +
  geom_area(position = "stack", colour = "black", size = 0.15, alpha = 0.85) +
  facet_wrap(~ core, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = module_colors_pal, name = "Module") +
  labs(
    title = "Module relative abundance through time",
    x = "Years Before Present, BP",
    y = "Relative abundance"
  ) +
  theme_ms +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 9, face = "bold")
  )


#########
# SUMMARY PLOTS
##########

## ── Plot B/C: module composition from current data ──────────────────────────

log_msg("Creating module composition plots...")

taxa_meta <- fread(file.path(RESULTS$stage1, "prokaryotes_taxa_metadata.tsv"))

taxa_mod <- merge(
  new_mods,
  taxa_meta[, .(taxon, functional_group, domain, phylum, signal_source)],
  by = "taxon",
  all.x = TRUE
)

# Preserve all module assignments; treat missing functional groups as Unknown
taxa_mod[is.na(functional_group), functional_group := "Unknown"]

# Order modules by taxon count for consistent plotting
module_counts <- taxa_mod[module != "grey", .N, by = module][order(-N)]
module_order <- module_counts$module

# Plot B: functional composition within each module
category_counts <- taxa_mod[module != "grey",
                            .N, by = .(module, functional_group)]
category_counts[, prop := N / sum(N), by = module]
category_counts[, module := factor(module, levels = module_order)]

group_order <- category_counts[, .(total = sum(N)), by = functional_group][order(-total)]
category_counts[, functional_group := factor(functional_group, levels = group_order$functional_group)]

p_module_function <- ggplot(category_counts, aes(x = module, y = prop, fill = functional_group)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_brewer(palette = "Set3", name = "Functional group") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.02)) +
  labs(
    title = "B. Functional composition by module",
    x = "Module",
    y = "Proportion of taxa"
  ) +
  theme_ms +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.key.size = unit(0.3, "cm")
  )

# Plot C: module size summary
module_counts <- module_counts[, module := factor(module, levels = module_order)]

module_colors <- c(
  turquoise = "#1B9E77",
  blue      = "#377EB8",
  brown     = "#A65628",
  yellow    = "#E6AB02",
  green     = "#66A61E",
  red       = "#E7298A",
  grey      = "#CCCCCC"
)

p_module_size <- ggplot(module_counts, aes(x = module, y = N, fill = module)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = N), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = module_colors) +
  labs(
    title = "C. Module sizes",
    x = "Module",
    y = "Number of taxa"
  ) +
  theme_ms +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )
