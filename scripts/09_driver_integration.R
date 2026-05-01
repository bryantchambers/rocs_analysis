#!/usr/bin/env Rscript
# 09_driver_integration.R — Identify Super-Driver taxa by combining Statistical and Topological importance
#
# Pipeline:
#   1. Load Fuzzy Forest importance and Network Statistics
#   2. Load Functional metrics (EMP thermodynamic capacity and TEA redox preference)
#   3. Join all datasets on taxon ID
#   4. Normalize importance scores to percentiles
#   5. Identify "Super-Drivers" (High Predictive Power AND High Topological Centrality)
#   6. Generate master integrated driver table and summary plots
#
# Outputs: results/importance/
#   integrated_driver_summary.tsv   Master table of all importance and functional metrics
#   super_driver_quadrant_plot.png  Plot comparing Statistical vs Topological importance

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

DRIVER_OUT <- RESULTS$importance
dir.create(DRIVER_OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading importance and network results...")
ff_imp <- fread("results/importance_fuzzy/ff_variable_importance.tsv")
net_stats <- fread("results/network_stats/network_metrics_summary.tsv")

log_msg("Loading functional metrics...")
emp_cap <- fread("results/emp/taxon_dg_capacity.tsv")
tea_pref <- fread("results/tea/oap_per_taxon.tsv")

# ── 2. Join Datasets ──────────────────────────────────────────────────────────

log_msg("Integrating results...")

# Start with network stats as it contains all taxa
integrated_dt <- copy(net_stats)

# Join Fuzzy Forest importance
# Note: FF only contains features that survived screening, so many will be NA
integrated_dt <- merge(integrated_dt, ff_imp[, .(taxon, variable_importance)], by = "taxon", all.x = TRUE)

# Join Functional metrics
integrated_dt <- merge(integrated_dt, emp_cap[, .(taxon, total_dg_capacity, primary_reaction)], by = "taxon", all.x = TRUE)
integrated_dt <- merge(integrated_dt, tea_pref[, .(taxon, oap_v3, dominant_class)], by = "taxon", all.x = TRUE)

# ── 3. Normalize & Rank ───────────────────────────────────────────────────────

log_msg("Calculating integrated importance ranks...")

# Handle NAs for ranking
integrated_dt[is.na(variable_importance), variable_importance := 0]

# Statistical Percentile
integrated_dt[, stat_percentile := rank(variable_importance) / .N]

# Composite Topological Percentile (Average of PageRank, Betweenness, and Vulnerability ranks)
integrated_dt[, topo_composite_percentile := (rank(pagerank) + rank(betweenness) + rank(vulnerability)) / (3 * .N)]

# Combined Score: Average of stat and topo percentiles
integrated_dt[, integrated_score := (stat_percentile + topo_composite_percentile) / 2]

# ── 4. Categorize Drivers (Tiered System) ─────────────────────────────────────

log_msg("Categorizing candidates into Tiers...")

integrated_dt[, driver_tier := fcase(
  # Tier 1: Top 10% in both
  stat_percentile >= 0.90 & topo_composite_percentile >= 0.90, "Tier 1: Super-Driver",
  
  # Tier 2: Top 20% in both OR specialized network role
  (stat_percentile >= 0.80 & topo_composite_percentile >= 0.80) | 
  is_potential_keystone == TRUE | is_hidden_gem == TRUE, "Tier 2: High Potential",
  
  # Tier 3: Statistical Specialists (Top 10% VIMP but lower topo)
  stat_percentile >= 0.90, "Tier 3: Predictive Specialist",
  
  default = "Peripheral"
)]

# Define Super-Drivers flag for plotting (Tier 1)
integrated_dt[, is_super_driver := (driver_tier == "Tier 1: Super-Driver")]

# ── 5. Save Results ───────────────────────────────────────────────────────────

log_msg(sprintf("Saving integrated driver table..."))
log_msg(sprintf("  Tier 1: %d", sum(integrated_dt$driver_tier == "Tier 1: Super-Driver")))
log_msg(sprintf("  Tier 2: %d", sum(integrated_dt$driver_tier == "Tier 2: High Potential")))
log_msg(sprintf("  Tier 3: %d", sum(integrated_dt$driver_tier == "Tier 3: Predictive Specialist")))

fwrite(integrated_dt, file.path(DRIVER_OUT, "integrated_driver_summary.tsv"), sep = "\t")

# ── 6. Visualization ──────────────────────────────────────────────────────────

log_msg("Generating comparison plots...")

p_quad <- ggplot(integrated_dt, aes(x = topo_composite_percentile, y = stat_percentile, color = driver_tier)) +
  geom_point(aes(size = integrated_score), alpha = 0.5) +
  scale_color_manual(values = c("Tier 1: Super-Driver" = "#D55E00", 
                               "Tier 2: High Potential" = "#0072B2", 
                               "Tier 3: Predictive Specialist" = "#009E73", 
                               "Peripheral" = "grey80")) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = "Statistical vs. Composite Topological Importance",
       subtitle = "Tiered Driver Identification System",
       x = "Composite Topological Percentile",
       y = "Statistical Percentile (Fuzzy Forest VIMP)") +
  geom_text_repel(data = integrated_dt[driver_tier != "Peripheral" & !is.na(species)], 
                  aes(label = species), size = 2.5, max.overlaps = 20)

ggsave(file.path(DRIVER_OUT, "super_driver_quadrant_plot.png"), p_quad, width = 12, height = 10)

log_msg("Done. Results in ", DRIVER_OUT)
