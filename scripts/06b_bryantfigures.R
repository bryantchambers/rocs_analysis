#!/usr/bin/env Rscript
# 06b_bryantfigures.R — evaluation and manuscript figure

#load packages and config
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
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

save_plots_to_results <- function(pattern = "^p_", suffix = "_draft", out_dir = FIGS) {
  plot_names <- ls(envir = .GlobalEnv, pattern = pattern)
  if (length(plot_names) == 0) {
    log_msg("No plots found matching pattern: ", pattern)
    return(NULL)
  }
  
  log_msg("Saving ", length(plot_names), " plots to ", out_dir)
  for (p_name in plot_names) {
    p <- get(p_name)
    if (inherits(p, "gg")) {
      file_base <- file.path(out_dir, paste0(p_name, suffix))
      ggsave(paste0(file_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)
      ggsave(paste0(file_base, ".pdf"), plot = p, width = 8, height = 6)
    }
  }
}

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


# ── LR04 δ¹⁸O reference panel ────────────────────────────────────────────────

log_msg("Creating LR04 climate reference panel...")

# Load MIS stage boundaries
mis_bounds <- fread(CLIMATE$mis_boundaries)
setnames(mis_bounds, c("stage", "start_kyr", "end_kyr", "climate"), 
         c("MIS", "start", "end", "climate"))
mis_bounds <- mis_bounds[end <= 150]

# Load LR04 benthic stack
lr04 <- fread(CLIMATE$lr04_stack)
setnames(lr04, c("Time_ka", "Benthic_d18O", "Standard_error"), 
         c("age_kyr", "d18O", "se"))
lr04_150 <- lr04[age_kyr <= 150]

# Label positions for MIS stages
mis_labels <- mis_bounds[start < 150]
mis_labels[, label_x := (start + pmin(end, 150)) / 2]
d18o_range <- range(lr04_150$d18O, na.rm = TRUE)

p_lr04_climate <- ggplot() +
  geom_rect(data = mis_bounds[start < 150],
            aes(xmin = start, xmax = pmin(end, 150), ymin = -Inf, ymax = Inf, fill = climate),
            alpha = 0.5) +
  scale_fill_manual(values = PALETTES$mis_climate, guide = "none") +
  geom_text(data = mis_labels,
            aes(x = label_x, y = d18o_range[1] - 0.15, label = MIS),
            size = 2.5, color = "grey30", fontface = "bold") +
  geom_line(data = lr04_150, aes(x = age_kyr, y = d18O), color = "black", linewidth = 0.5) +
  scale_x_reverse(limits = c(150, 0), breaks = seq(150, 0, -25), expand = c(0.01, 0)) +
  scale_y_reverse() +
  labs(x = NULL, y = expression(delta^{18}*"O (‰)"), 
       title = "LR04 benthic δ¹⁸O stack") +
  theme_ms +
  theme(axis.text.x = element_blank(), 
        plot.margin = margin(5, 5, 0, 5),
        axis.title.y = element_text(size = 9))

log_msg("LR04 panel created")

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
  facet_wrap(~ core, scales = "fixed", ncol = 1) +
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
log_msg("Module abundance plot created")
plot(p_module_abund)

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

log_msg("Module composition plots created") 



# =========================================
# Plot D: state activity through time
# =========================================

# ── HMM states area plot ─────────────────────────────────────────────────────

log_msg("Creating HMM states plot...")

hmm_states <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))

# Merge to get y_bp from metadata
hmm_states <- merge(hmm_states, meta[, .(label, y_bp)], by.x = "sample", by.y = "label")
hmm_states[, age_kyr := y_bp / 1000]

# Order by core and y_bp
hmm_states <- hmm_states[order(core, y_bp)]

# Define colors for states (similar to module colors)
state_colors <- c(
  "1" = "#D51414", # Turquoise
  "2" = "#377EB8", # Blue
  "3" = "#A65628", # Brown
  "4" = "#E6AB02", # Yellow
  "5" = "#1e20a6"  # Green
)

p_hmm_states <- ggplot(hmm_states, aes(x = age_kyr, y = 1, fill = factor(state))) +
  # Width is set large (15 ka) to ensure tiles overlap and form a continuous bar 
  # despite irregular sampling intervals.
  geom_tile(width = 15, height = 1, colour = NA) +
  scale_fill_manual(values = state_colors, name = "State") +
  scale_x_continuous(limits = c(0, 600)) +
  facet_wrap(~ core, ncol = 1) +
  labs(
    title = "HMM States through time",
    x = "Age (ka)",
    y = NULL
  ) +
  theme_ms +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(0.2, "lines")
  )

log_msg("HMM states plot created")
plot(p_hmm_states)

# ── Save all figures ──────────────────────────────────────────────────────────

log_msg("Exporting figures...")


save_plots_to_results()

