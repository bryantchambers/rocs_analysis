#!/usr/bin/env Rscript
# 14_driver_quadrants.R — Visualize Statistical vs. Topological importance
#
# Pipeline:
#   1. Load Integrated Driver Summary
#   2. Generate polished quadrant plot with Tier highlights
#
# Outputs: results/figures/
#   driver_quadrant_analysis.png

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

FIG_OUT <- RESULTS$figures

log_msg("Loading master driver data...")
master <- fread(file.path(RESULTS$importance, "functional_driver_master.tsv"))

p_quad <- ggplot(master, aes(x = topo_composite_percentile, y = stat_percentile, color = driver_tier)) +
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
       x = "Composite Topological Influence (PageRank + Vulnerability + Betweenness)",
       y = "Statistical Predictive Importance (Fuzzy Forest)",
       color = "Driver Tier",
       size = "Integrated Score") +
  theme(legend.position = "right") +
  geom_text_repel(data = master[driver_tier != "Peripheral" & (integrated_score > 0.8 | is_super_driver == TRUE)], 
                  aes(label = species), size = 3, max.overlaps = 30)

ggsave(file.path(FIG_OUT, "driver_quadrant_analysis.png"), p_quad, width = 12, height = 10)

log_msg("Done. Results in ", FIG_OUT)
