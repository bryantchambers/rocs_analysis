#!/usr/bin/env Rscript
# 13_state_transition_network.R — Visualize module-level connectivity shifts between states
#
# Pipeline:
#   1. Load state-specific connectivity stats and inter-module edges
#   2. Construct a module-level "Meta-Network" for each HMM state
#   3. Edge weight = Density of connections between modules
#   4. Node size = Module size or average state-specific abundance
#   5. Generate a multi-panel network plot (Glacial vs Interglacial)
#
# Outputs: results/figures/
#   state_transition_meta_network.png    The "storyboard" of module interactions

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(WGCNA)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

FIG_OUT <- RESULTS$figures

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading connectivity and module data...")
vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
mods <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
states <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))

common_samples <- intersect(rownames(vst), states$sample)
vst_sub <- vst[common_samples, ]
states_sub <- states[match(common_samples, sample), ]

soft_power <- 12
adj_global <- adjacency(vst_sub, power = soft_power, type = "signed")
adj_global[adj_global < 0.1] <- 0

# ── 2. Calculate Module-Level Edges per State ─────────────────────────────────

unique_states <- unique(states_sub$label)
mod_colors <- c("turquoise" = "turquoise", "blue" = "blue", "brown" = "brown", 
                "yellow" = "yellow", "green" = "green", "red" = "red")
all_mods <- names(mod_colors)

# We'll use a fixed layout for all states
g_layout_base <- graph_from_literal(turquoise-blue-brown-yellow-green-red-turquoise)
layout_mat <- layout_in_circle(g_layout_base)
rownames(layout_mat) <- V(g_layout_base)$name

png(file.path(FIG_OUT, "state_transition_meta_network.png"), width = 1600, height = 1200, res = 150)
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))

for (st in unique_states) {
  log_msg(sprintf("  Drawing meta-network for: %s", st))
  
  st_samples <- states_sub[label == st, sample]
  st_vst <- vst_sub[st_samples, ]
  avg_abund <- colMeans(st_vst)
  active_taxa <- names(avg_abund)[avg_abund > quantile(avg_abund, 0.25)]
  
  st_adj <- adj_global[active_taxa, active_taxa]
  g_taxa <- graph_from_adjacency_matrix(st_adj, weighted=TRUE, mode="undirected", diag=FALSE)
  
  # Contract vertices by module
  V(g_taxa)$module <- mods[match(V(g_taxa)$name, taxon), module]
  g_taxa <- delete_vertices(g_taxa, V(g_taxa)[is.na(module) | module == "grey"])
  
  # Manual contraction to avoid igraph vertex mapping issues
  m_adj <- matrix(0, nrow=length(all_mods), ncol=length(all_mods), dimnames=list(all_mods, all_mods))
  edges_dt <- as.data.table(as_data_frame(g_taxa))
  edges_dt[, from_mod := mods[match(from, taxon), module]]
  edges_dt[, to_mod := mods[match(to, taxon), module]]
  
  meta_summary <- edges_dt[, .(w = sum(weight)), by = .(from_mod, to_mod)]
  for(i in 1:nrow(meta_summary)) {
    m_adj[meta_summary$from_mod[i], meta_summary$to_mod[i]] <- meta_summary$w[i]
  }
  
  g_meta <- graph_from_adjacency_matrix(m_adj, weighted=TRUE, mode="undirected", diag=FALSE)
  
  # Plotting
  V(g_meta)$color <- mod_colors[V(g_meta)$name]
  E(g_meta)$width <- log10(E(g_meta)$weight + 1) * 2
  
  plot(g_meta, 
       layout = layout_mat[V(g_meta)$name, ], 
       vertex.size = 40,
       vertex.label.color = "black",
       vertex.label.cex = 1.2,
       vertex.label.dist = 0,
       edge.color = "grey60",
       edge.alpha = 0.5,
       main = paste("State:", st))
}
dev.off()

log_msg("Done. Results in ", FIG_OUT)

log_msg("Done. Results in ", FIG_OUT)
