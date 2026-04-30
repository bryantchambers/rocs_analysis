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

# Percentiles (0 to 1 scale, higher is more important)
integrated_dt[, stat_percentile := rank(variable_importance) / .N]
integrated_dt[, topo_percentile := rank(pagerank) / .N] # Using PageRank as a proxy for global topological influence

# Combined Score: Average of percentiles
integrated_dt[, integrated_score := (stat_percentile + topo_percentile) / 2]

# Define Super-Drivers: Top 5% in BOTH or Top 2.5% combined
integrated_dt[, is_super_driver := (stat_percentile >= 0.95 & topo_percentile >= 0.95) | (integrated_score >= 0.975)]

# ── 4. Categorize Drivers ─────────────────────────────────────────────────────

integrated_dt[, driver_category := fcase(
  stat_percentile >= 0.90 & topo_percentile >= 0.90, "Super-Driver",
  stat_percentile >= 0.90 & topo_percentile < 0.90,  "Predictive Specialist",
  stat_percentile < 0.90 & topo_percentile >= 0.90,  "Topological Keystone",
  default = "Peripheral"
)]

# ── 5. Save Results ───────────────────────────────────────────────────────────

log_msg(sprintf("Saving integrated driver table (%d super-drivers identified)...", sum(integrated_dt$is_super_driver)))
fwrite(integrated_dt, file.path(DRIVER_OUT, "integrated_driver_summary.tsv"), sep = "\t")

# ── 6. Visualization ──────────────────────────────────────────────────────────

log_msg("Generating comparison plots...")

p_quad <- ggplot(integrated_dt, aes(x = topo_percentile, y = stat_percentile, color = module)) +
  geom_point(aes(size = integrated_score), alpha = 0.4) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = "Statistical vs. Topological Driver Importance",
       subtitle = "Top-right quadrant indicates Super-Drivers",
       x = "Topological Percentile (PageRank)",
       y = "Statistical Percentile (Fuzzy Forest VIMP)") +
  geom_text_repel(data = integrated_dt[is_super_driver == TRUE & !is.na(species)], 
                  aes(label = species), size = 3, max.overlaps = 15)

ggsave(file.path(DRIVER_OUT, "super_driver_quadrant_plot.png"), p_quad, width = 12, height = 10)

log_msg("Done. Results in ", DRIVER_OUT)
