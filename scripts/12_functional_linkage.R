#!/usr/bin/env Rscript
# 12_functional_linkage.R — Link drivers and bridge taxa to metabolic/redox traits
#
# Pipeline:
#   1. Load Integrated Drivers, Climate Sensitivity, and State Bridges
#   2. Load Functional Annotations (EMP thermodynamic capacity and TEA redox classes)
#   3. Calculate functional enrichment for each Driver Tier and HMM State Bridge set
#   4. Identify "Functional Hubs": High connectivity + High metabolic capacity
#   5. Map "Redox Bridges": Taxa bridging modules while carrying specific TEA pathways
#
# Outputs: results/importance/
#   functional_driver_master.tsv      Drivers annotated with all functional and climate traits
#   state_functional_enrichment.tsv   Enrichment of TEA classes in state-specific bridge sets
#   functional_linkage_summary.png    Plot showing EMP vs. Connectivity colored by TEA

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

OUT_DIR <- RESULTS$importance

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading integrated results and functional traits...")
drivers <- fread(file.path(RESULTS$importance, "integrated_driver_summary.tsv"))
climate <- fread(file.path(RESULTS$importance, "climate_sensitivity_results.tsv"))
bridges <- fread(file.path(RESULTS$network_stats, "bridge_taxa_by_state.tsv"))

# Functional traits were already joined in 09, but let's make sure we have the full set
# integrated_driver_summary.tsv contains: total_dg_capacity, primary_reaction, oap_v3, dominant_class

# ── 2. Functional Profile of Tiers ────────────────────────────────────────────

log_msg("Summarizing functional profiles by Driver Tier...")

# Calculate % of each Tier carrying specific TEA classes
tier_tea_summary <- drivers[driver_tier != "Peripheral", .(
  n_taxa = .N,
  avg_emp = mean(total_dg_capacity, na.rm = TRUE),
  prop_aerobic = sum(dominant_class == "O2", na.rm = TRUE) / .N,
  prop_sulfate = sum(dominant_class == "SO4", na.rm = TRUE) / .N,
  prop_methane = sum(dominant_class == "CH4", na.rm = TRUE) / .N
), by = driver_tier]

# ── 3. State-Specific Bridge Functionality ───────────────────────────────────

log_msg("Analyzing functional roles of State-Specific Bridges...")

# Join bridges with functional data
bridge_func <- merge(bridges, drivers[, .(taxon, total_dg_capacity, primary_reaction, oap_v3, dominant_class)], by = "taxon")

state_func_enrich <- bridge_func[, .(
  n_bridges = .N,
  avg_emp = mean(total_dg_capacity, na.rm = TRUE),
  top_tea = names(sort(table(dominant_class), decreasing = TRUE))[1]
), by = state]

# ── 4. Identify Functional Super-Hubs ────────────────────────────────────────

log_msg("Identifying Functional Super-Hubs...")

# Relaxing criteria to ensure we capture the intersection of connectivity and capacity
# Functional Super-Hub = Top 25% Integrated Score AND Above-median EMP capacity
score_thresh <- quantile(drivers$integrated_score, 0.75, na.rm = TRUE)
emp_thresh <- median(drivers$total_dg_capacity, na.rm = TRUE)

drivers[, is_functional_hub := (integrated_score >= score_thresh & total_dg_capacity >= emp_thresh)]

log_msg(sprintf("  Identified %d Functional Super-Hubs using score >= %.2f and EMP >= %.2f", 
                sum(drivers$is_functional_hub, na.rm = TRUE), score_thresh, emp_thresh))

# ── 5. Save Integrated Master Table ──────────────────────────────────────────

# Merge climate sensitivity into the master driver table
master_drivers <- merge(drivers, climate[, .(taxon, coef_d18O, fdr_d18O)], by = "taxon", all.x = TRUE)

# Handle cases where drivers might be missing species names
master_drivers[is.na(species) | species == "", species := taxon]

fwrite(master_drivers, file.path(OUT_DIR, "functional_driver_master.tsv"), sep = "\t")
fwrite(state_func_enrich, file.path(OUT_DIR, "state_functional_enrichment.tsv"), sep = "\t")

# ── 6. Visualization ──────────────────────────────────────────────────────────

log_msg("Generating functional linkage plots...")

# Plot 1: Topological Influence vs. Metabolic Capacity (EMP)
p_func <- ggplot(master_drivers[driver_tier != "Peripheral" & !is.na(total_dg_capacity)], 
                aes(x = topo_composite_percentile, y = total_dg_capacity, color = dominant_class)) +
  geom_point(aes(size = integrated_score), alpha = 0.7) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1", na.value = "grey80") +
  labs(title = "Topological Influence vs. Metabolic Capacity",
       subtitle = "Points represent top candidate drivers; Size = Integrated Importance",
       x = "Composite Topological Percentile",
       y = "Total DG Capacity (EMP)",
       color = "TEA Class") +
  geom_text_repel(data = master_drivers[is_functional_hub == TRUE & driver_tier != "Peripheral"], 
                  aes(label = species), size = 2.5, max.overlaps = 20)

ggsave(file.path(OUT_DIR, "functional_linkage_connectivity.png"), p_func, width = 12, height = 10)

# Plot 2: Climate Sensitivity (Glacial vs Interglacial) vs. Metabolic Capacity
p_clim_func <- ggplot(master_drivers[!is.na(coef_d18O) & driver_tier != "Peripheral" & !is.na(total_dg_capacity)], 
                     aes(x = coef_d18O, y = total_dg_capacity, color = dominant_class)) +
  geom_point(aes(size = -log10(fdr_d18O)), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1", na.value = "grey80") +
  labs(title = "Climate Sensitivity vs. Metabolic Capacity",
       subtitle = "Positive Coef = Glacial Preference; Size = Significance",
       x = "GLS Coefficient (d18O Influence)",
       y = "Total DG Capacity (EMP)",
       color = "TEA Class") +
  geom_text_repel(data = master_drivers[abs(coef_d18O) > 0.8 & driver_tier != "Peripheral"], 
                  aes(label = species), size = 2.5, max.overlaps = 20)

ggsave(file.path(OUT_DIR, "functional_climate_sensitivity.png"), p_clim_func, width = 12, height = 10)

log_msg("Done. Results in ", OUT_DIR)
