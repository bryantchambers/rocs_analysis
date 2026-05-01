#!/usr/bin/env Rscript
# 16_final_story_visualization.R — Integrated Timeline of States, Function, and Drivers
#
# Pipeline:
#   1. Load HMM States, Climate (LR04), Functional Indices (EMP/TEA), and Bridge Taxa
#   2. Panel A: Climate (d18O) + HMM States Timeline
#   3. Panel B: Functional Index shifts (EMP and OAP) along timeline
#   4. Panel C: "Bridge Intensity" (Inter-module connectivity) along timeline
#   5. Panel D: Functional breakdown of top bridges at transition points
#
# Outputs: results/figures/
#   final_biological_narrative_timeline.png

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

FIG_OUT <- RESULTS$figures

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading integrated datasets...")

# Climate (LR04)
# Using fill=TRUE to handle potential trailing whitespace/formatting issues noted in earlier run
lr04 <- fread("old/Lisiecki2005_copy.txt", fill = TRUE)
# The file has 3 columns: Time_ka, Benthic_d18O, Standard_error
setnames(lr04, 1:2, c("age_kyr", "d18O")) 
lr04 <- lr04[, .(age_kyr, d18O)]

# HMM States
states <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))

# Functional Indices (Per Sample)
emp_sample <- fread(file.path(RESULTS$emp, "emp_sap_per_sample.tsv"))
tea_sample <- fread(file.path(RESULTS$tea, "tea_indices_per_sample.tsv"))

# Bridge Taxa & Master Drivers
bridges <- fread(file.path(RESULTS$network_stats, "bridge_taxa_by_state.tsv"))
master <- fread(file.path(RESULTS$importance, "functional_driver_master.tsv"))

# Sample Metadata
meta <- fread(UPSTREAM$metadata)
meta_clean <- meta[, .(sample = label, core, age_kyr = y_bp / 1000)]

# ── 2. Prepare Combined Timeline Data ─────────────────────────────────────────

log_msg("Merging timeline data...")

# Merge everything on sample/age
timeline_dt <- merge(meta_clean, states[, .(sample, state_label = label)], by = "sample")
timeline_dt <- merge(timeline_dt, emp_sample[, .(sample, EMP)], by = "sample", all.x = TRUE)
timeline_dt <- merge(timeline_dt, tea_sample[, .(sample, OAP)], by = "sample", all.x = TRUE)

# Calculate Module Relative Abundances (Logic from 06b)
log_msg("Calculating module relative abundances for timeline...")
prok_raw <- fread(UPSTREAM$tax_damage)
prok <- prok_raw[is_dmg == "Damaged" & label %in% meta_clean$sample]
prok_agg <- prok[, .(n_reads = sum(n_reads)), by = .(subspecies, label)]
mods_map <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
setnames(mods_map, "taxon", "subspecies")
prok_mod <- merge(prok_agg, mods_map, by = "subspecies", all.x = TRUE)
prok_mod[is.na(module), module := "grey"]
prok_mod[, total_reads := sum(n_reads), by = label]
prok_mod[, rel_abund := n_reads / total_reads]

module_abund <- prok_mod[, .(rel_abund = sum(rel_abund)), by = .(module, sample = label)]
module_abund <- merge(module_abund, meta_clean[, .(sample, age_kyr)], by = "sample")

# Filter for the last 150kyr
timeline_sub <- timeline_dt[age_kyr <= 150]
lr04_sub <- lr04[age_kyr <= 150]
module_sub <- module_abund[age_kyr <= 150]

# ── 3. Panel A: Climate & States ──────────────────────────────────────────────

log_msg("Creating Panel A: Climate & States...")

state_colors <- c("G-A" = "#D55E00", "IG-B" = "#0072B2", "IG-C" = "#009E73", 
                 "IG-E" = "#CC79A7", "IG-A" = "#F0E442")

timeline_sub <- timeline_sub[order(age_kyr)]
timeline_sub[, next_age := shift(age_kyr, type = "lead")]
timeline_sub[is.na(next_age), next_age := age_kyr + 1]

p_clim <- ggplot() +
  geom_line(data = lr04_sub, aes(x = age_kyr, y = d18O), color = "grey70", alpha = 0.8) +
  geom_rect(data = timeline_sub, 
            aes(xmin = age_kyr, xmax = next_age, ymin = 5.2, ymax = 5.5, fill = state_label),
            color = NA) +
  scale_fill_manual(values = state_colors) +
  scale_y_reverse(limits = c(5.5, 3.0)) + 
  scale_x_reverse(limits = c(150, 0), breaks = seq(0, 150, 25)) +
  theme_minimal() +
  labs(y = expression(delta^{18}*"O (‰)"), x = NULL, 
       title = "A. Climate Forcing & Ecological States", fill = "HMM State") +
  theme(legend.position = "right", panel.grid.minor = element_blank(), axis.text.x = element_blank())

# ── 4. Panel B: Module Relative Abundance ─────────────────────────────────────

log_msg("Creating Panel B: Module Abundances...")

module_colors_pal <- c(
  turquoise = "#1B9E77", blue = "#377EB8", brown = "#A65628",
  yellow = "#E6AB02", green = "#66A61E", red = "#E7298A", grey = "#CCCCCC"
)

# Smoothing module abundances for a cleaner area plot
p_mod_abund <- ggplot(module_sub[module != "grey"], aes(x = age_kyr, y = rel_abund, fill = module)) +
  geom_area(position = "stack", alpha = 0.8, color = "white", size = 0.1) +
  scale_fill_manual(values = module_colors_pal) +
  scale_x_reverse(limits = c(150, 0), breaks = seq(0, 150, 25)) +
  theme_minimal() +
  labs(y = "Rel. Abund.", x = NULL, title = "B. Community Composition (WGCNA Modules)", fill = "Module") +
  theme(legend.position = "right", axis.text.x = element_blank())

# ── 5. Panel C: Functional Shifts (EMP & OAP) ─────────────────────────────────

log_msg("Creating Panel C: Functional Shifts...")

p_func <- ggplot(timeline_sub, aes(x = age_kyr)) +
  geom_smooth(aes(y = EMP / max(EMP, na.rm=T), color = "EMP (Metabolic Capacity)"), method = "loess", se = FALSE) +
  geom_smooth(aes(y = (OAP - min(OAP, na.rm=T)) / (max(OAP, na.rm=T) - min(OAP, na.rm=T)), 
                  color = "OAP (Redox Potential)"), method = "loess", se = FALSE) +
  scale_x_reverse(limits = c(150, 0), breaks = seq(0, 150, 25)) +
  theme_minimal() +
  labs(y = "Norm. Index", x = NULL, title = "C. Functional Response (Metabolic & Redox)", color = "Index") +
  theme(legend.position = "right", axis.text.x = element_blank())

# ── 6. Panel D: Network Coordination (Bridge Intensity) ───────────────────────

log_msg("Creating Panel D: Network Coordination...")

sulfur_taxa <- master[dominant_class == "SO4", taxon]
vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
sulfur_abund <- data.table(sample = rownames(vst), 
                          sulfur_intensity = rowSums(vst[, intersect(colnames(vst), sulfur_taxa), drop=FALSE]))

timeline_sub <- merge(timeline_sub, sulfur_abund, by = "sample", all.x = TRUE)

p_bridge <- ggplot(timeline_sub, aes(x = age_kyr, y = sulfur_intensity)) +
  geom_smooth(method = "loess", color = "#CC79A7", fill = "#CC79A7", alpha = 0.2) +
  geom_point(alpha = 0.3) +
  scale_x_reverse(limits = c(150, 0), breaks = seq(0, 150, 25)) +
  theme_minimal() +
  labs(y = "CLR Abundance", x = "Age (ka)", 
       title = "D. Redox Bridge Activity (Sulfur Metabolism)",
       subtitle = "Aggregate abundance of top sulfur-reducing bridge taxa") +
  theme(panel.grid.minor = element_blank())

# ── 7. Panel E: Key Drivers & Metabolic Story ─────────────────────────────────

log_msg("Creating Panel E: Top Drivers Table...")

top_story_taxa <- master[is_functional_hub == TRUE | driver_tier == "Tier 1: Super-Driver"]
top_story_taxa <- top_story_taxa[order(-integrated_score)]
story_table <- top_story_taxa[1:10, .(Species = species, 
                                     Metabolism = functional_group,
                                     Tier = driver_tier,
                                     Climate = ifelse(coef_d18O > 0, "Glacial", "Interglacial"))]

p_table <- tableGrob(story_table, rows = NULL, theme = ttheme_minimal(base_size = 7))

# ── 8. Assemble and Save ──────────────────────────────────────────────────────

log_msg("Assembling final multi-panel figure...")

final_p <- (p_clim / p_mod_abund / p_func / p_bridge / wrap_elements(p_table)) + 
  plot_layout(heights = c(1, 1, 1, 1, 0.8)) +
  plot_annotation(title = "The Biological Engine of the Ancient Ocean",
                  subtitle = "Integrated view of climate cycles, ecological states, community shifts, and functional drivers",
                  caption = "Phase 16: Final Integrated Narrative")

ggsave(file.path(FIG_OUT, "final_biological_narrative_timeline.png"), final_p, width = 12, height = 18)

log_msg("Done. Results in ", FIG_OUT)
