#!/usr/bin/env Rscript
# 11_state_networks.R — Extract and analyze state-specific sub-networks
#
# Pipeline:
#   1. Load Integrated Drivers and Climate Sensitivity
#   2. Load TOM and HMM State assignments
#   3. For each HMM state:
#      a. Identify "Active" taxa (mean abundance > 0)
#      b. Extract sub-network from TOM
#      c. Calculate inter-module connectivity
#      d. Identify "Bridge Taxa" (Nodes with high edges to other modules)
#   4. Compare connectivity patterns across Glacial vs Interglacial states
#
# Outputs: results/network_stats/
#   state_connectivity_summary.tsv    Inter-module edge density per state
#   bridge_taxa_by_state.tsv          Top bridging taxa for each state
#   state_network_stats.tsv           Global network metrics (density, transitivity) per state

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(igraph)
  library(WGCNA)
  library(ggplot2)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

NET_OUT <- RESULTS$network_stats

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading TOM and driver metadata...")
# We'll re-calculate a slightly more filtered TOM if needed, or use the global one
vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
mods <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
states <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))
drivers <- fread(file.path(RESULTS$importance, "integrated_driver_summary.tsv"))

# Ensure sample alignment
common_samples <- intersect(rownames(vst), states$sample)
vst_sub <- vst[common_samples, ]
states_sub <- states[match(common_samples, sample), ]

# Re-calculate adjacency (global) to get full weights
soft_power <- 12 # From script 08
adj_global <- adjacency(vst_sub, power = soft_power, type = "signed")

# ── 2. State-Specific Analysis ────────────────────────────────────────────────

unique_states <- unique(states_sub$label)
log_msg(sprintf("Analyzing sub-networks for %d HMM states: %s", 
                length(unique_states), paste(unique_states, collapse=", ")))

state_metrics <- list()
bridge_list <- list()

for (st in unique_states) {
  log_msg(sprintf("  Processing state: %s", st))
  
  # Identify samples in this state
  st_samples <- states_sub[label == st, sample]
  st_vst <- vst_sub[st_samples, ]
  
  # Active taxa: defined as having abundance in this state (or top candidates)
  # For structural comparison, we use the same node set but look at edge shifts?
  # Actually, WGCNA TOM is static. To see "state-specific" networks, 
  # we usually filter for taxa present/abundant in that state.
  
  avg_abund <- colMeans(st_vst)
  active_taxa <- names(avg_abund)[avg_abund > quantile(avg_abund, 0.25)] # Keep top 75% active
  
  # Extract sub-adjacency
  sub_adj <- adj_global[active_taxa, active_taxa]
  
  # Construct graph
  sub_adj[sub_adj < 0.1] <- 0 # Sparsify for stats
  g_st <- graph_from_adjacency_matrix(sub_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # ── 3. Calculate Bridging/Inter-module ──────────────────────────────────────
  
  # Mapping nodes to modules
  node_mods <- mods[match(V(g_st)$name, taxon), module]
  
  # Inter-module connectivity: How many edges connect different modules?
  edges_dt <- as.data.table(as_data_frame(g_st, what = "edges"))
  edges_dt[, from_mod := node_mods[match(from, V(g_st)$name)]]
  edges_dt[, to_mod := node_mods[match(to, V(g_st)$name)]]
  edges_dt[, is_inter_module := (from_mod != to_mod)]
  
  inter_mod_ratio <- sum(edges_dt$is_inter_module) / nrow(edges_dt)
  
  # Identify "Bridge Taxa": Nodes with the most inter-module edges
  bridge_taxa <- edges_dt[is_inter_module == TRUE, .(inter_edges = .N, inter_weight = sum(weight)), by = .(taxon = from)]
  # Add the 'to' side as well
  bridge_taxa_to <- edges_dt[is_inter_module == TRUE, .(inter_edges = .N, inter_weight = sum(weight)), by = .(taxon = to)]
  bridge_combined <- rbind(bridge_taxa, bridge_taxa_to)[, .(inter_edges = sum(inter_edges), inter_weight = sum(inter_weight)), by = taxon]
  
  setorder(bridge_combined, -inter_weight)
  bridge_combined[, state := st]
  bridge_list[[st]] <- head(bridge_combined, 50) # Top 50 bridges per state
  
  # Global state stats
  state_metrics[[st]] <- data.table(
    state = st,
    n_nodes = vcount(g_st),
    n_edges = ecount(g_st),
    density = edge_density(g_st),
    transitivity = transitivity(g_st),
    inter_mod_ratio = inter_mod_ratio
  )
}

# ── 4. Save and Summarize ─────────────────────────────────────────────────────

log_msg("Saving state-specific results...")
final_metrics <- rbindlist(state_metrics)
final_bridges <- rbindlist(bridge_list)

# Annotate bridges with species and driver tier
final_bridges <- merge(final_bridges, drivers[, .(taxon, species, module, driver_tier)], by = "taxon")
setorder(final_bridges, state, -inter_weight)

fwrite(final_metrics, file.path(NET_OUT, "state_network_stats.tsv"), sep = "\t")
fwrite(final_bridges, file.path(NET_OUT, "bridge_taxa_by_state.tsv"), sep = "\t")

# ── 5. Visualization ──────────────────────────────────────────────────────────

log_msg("Generating comparison plots...")

p_density <- ggplot(final_metrics, aes(x = state, y = inter_mod_ratio, fill = state)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Inter-Module Connectivity by HMM State",
       y = "Ratio of Inter-Module Edges")

ggsave(file.path(NET_OUT, "state_connectivity_comparison.png"), p_density, width = 8, height = 6)

log_msg("Done. Results in ", NET_OUT)
