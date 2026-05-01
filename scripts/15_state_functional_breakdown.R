#!/usr/bin/env Rscript
# 15_state_functional_breakdown.R — Visualize functional composition of bridge taxa per state
#
# Pipeline:
#   1. Load State Bridge Taxa and Functional traits
#   2. Generate stacked bar plots of TEA classes per state
#
# Outputs: results/figures/
#   state_bridge_functional_composition.png

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

FIG_OUT <- RESULTS$figures

log_msg("Loading bridge and functional data...")
bridges <- fread(file.path(RESULTS$network_stats, "bridge_taxa_by_state.tsv"))
drivers <- fread(file.path(RESULTS$importance, "functional_driver_master.tsv"))

# Merge functionality into bridges
bridge_func <- merge(bridges, drivers[, .(taxon, dominant_class)], by = "taxon")

# Handle NAs
bridge_func[is.na(dominant_class), dominant_class := "Unknown"]

# Count occurrences of each TEA class per state
composition <- bridge_func[, .(count = .N), by = .(state, dominant_class)]

p_comp <- ggplot(composition, aes(x = state, y = count, fill = dominant_class)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Functional Composition of State Bridge Taxa",
       subtitle = "Top 50 bridging taxa per HMM state",
       x = "HMM State",
       y = "Proportion of Bridge Taxa",
       fill = "TEA Class (Redox)") +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_OUT, "state_bridge_functional_composition.png"), p_comp, width = 10, height = 8)

log_msg("Done. Results in ", FIG_OUT)
