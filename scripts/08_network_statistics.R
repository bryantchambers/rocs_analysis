#!/usr/bin/env Rscript
# 08_network_statistics.R — Compute network topology and centrality metrics
#
# Pipeline:
#   1. Load WGCNA consensus network and module assignments
#   2. Construct igraph object from Topological Overlap Matrix (TOM)
#   3. Calculate Within-Module Degree (z) and Participation Coefficient (p)
#   4. Calculate Centrality metrics (PageRank, Closeness, Betweenness)
#   5. Calculate Bridging Centrality and Vulnerability (Nodal Efficiency drop)
#   6. Annotate with taxonomic metadata and save results
#
# Outputs: results/network_stats/
#   network_metrics_summary.tsv   Master table of all nodal metrics
#   zp_plot_all.png               Z-P plot for the full network
#   zp_plot_by_module.png         Z-P plots faceted by module

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
  library(igraph)
  library(ggplot2)
  library(ggrepel)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

STATS_OUT <- RESULTS$network_stats
dir.create(STATS_OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading WGCNA results...")
net <- readRDS("results/stage1/wgcna/consensus_wgcna.rds")
mods <- fread("results/stage1/wgcna/module_assignments.tsv")
tax_meta <- fread("results/stage1/prokaryotes_taxa_metadata.tsv")

# Get proper species/genus names from raw data
log_msg("Extracting taxon names from raw damage summary...")
tax_map <- unique(fread(UPSTREAM$tax_damage, select = c("subspecies", "genus", "species")))
setnames(tax_map, "subspecies", "taxon")
# Keep only one name per taxon if duplicates exist
tax_map <- tax_map[, .(genus = genus[1], species = species[1]), by = taxon]

# Re-loading VST data to calculate TOM
vst <- readRDS("results/stage1/prokaryotes_vst.rds")
common_taxa <- mods$taxon
vst_sub <- vst[, common_taxa]

log_msg(sprintf("Processing %d taxa...", length(common_taxa)))

# ── 2. Construct Network ──────────────────────────────────────────────────────

log_msg("Calculating Adjacency and TOM...")
soft_power <- net$power
if(is.null(soft_power)) soft_power <- 12 

# Calculate adjacency (signed)
adj <- adjacency(vst_sub, power = soft_power, type = "signed")

# Calculate TOM
tom <- TOMsimilarity(adj, TOMType = "signed")
rownames(tom) <- colnames(tom) <- common_taxa

# Convert to igraph
g <- graph_from_adjacency_matrix(tom, mode = "undirected", weighted = TRUE, diag = FALSE)

# ── 3. Z-P Metrics ────────────────────────────────────────────────────────────

log_msg("Calculating Z and P scores...")

calc_z_p <- function(adj_mat, module_assigns) {
  taxa <- rownames(adj_mat)
  unique_mods <- unique(module_assigns$module)
  
  k_total <- rowSums(adj_mat)
  k_within <- matrix(0, nrow = length(taxa), ncol = length(unique_mods))
  colnames(k_within) <- unique_mods
  rownames(k_within) <- taxa
  
  for (m in unique_mods) {
    m_taxa <- module_assigns[module == m, taxon]
    if (length(m_taxa) > 1) {
      k_within[, m] <- rowSums(adj_mat[, m_taxa])
    } else if (length(m_taxa) == 1) {
      k_within[, m] <- adj_mat[, m_taxa]
    }
  }
  
  ki_mi <- sapply(taxa, function(x) {
    m <- module_assigns[taxon == x, module]
    k_within[x, m]
  })
  
  z_scores <- sapply(taxa, function(x) {
    m <- module_assigns[taxon == x, module]
    m_taxa <- module_assigns[module == m, taxon]
    m_ks <- ki_mi[m_taxa]
    if (length(m_ks) < 2 || sd(m_ks) == 0) return(0)
    (ki_mi[x] - mean(m_ks)) / sd(m_ks)
  })
  
  p_scores <- 1 - rowSums((k_within / k_total)^2)
  return(data.table(taxon = taxa, z_score = z_scores, p_score = p_scores))
}

zp_dt <- calc_z_p(adj, mods)

zp_dt[, role := fcase(
  z_score >= 2.5 & p_score <= 0.62, "Module Hub",
  z_score >= 2.5 & p_score > 0.62,  "Network Hub",
  z_score < 2.5  & p_score > 0.62,  "Connector",
  z_score < 2.5  & p_score <= 0.62, "Peripheral"
)]

# ── 4. Centrality Metrics ─────────────────────────────────────────────────────

log_msg("Calculating Centrality (PageRank, Closeness, Betweenness)...")

# For distance-based metrics, we use 1 - TOM or 1/TOM as distance
# Betweenness and Closeness use weights as lengths (shorter is better)
dist_weights <- 1 - E(g)$weight
closeness_val <- closeness(g, weights = dist_weights, normalized = TRUE)
betweenness_val <- betweenness(g, weights = dist_weights, normalized = TRUE)
pagerank_val <- page_rank(g, weights = E(g)$weight)$vector # PageRank uses weights as strength

centrality_dt <- data.table(
  taxon = names(closeness_val),
  closeness = closeness_val,
  betweenness = betweenness_val,
  pagerank = pagerank_val
)

# ── 5. Bridging & Vulnerability ───────────────────────────────────────────────

log_msg("Calculating Bridging Centrality...")

calc_bridging_coeff <- function(graph) {
  deg <- degree(graph)
  bc <- sapply(V(graph), function(v) {
    neighs <- neighbors(graph, v)
    if (length(neighs) == 0) return(0)
    (1/deg[v]) / sum(1/deg[neighs])
  })
  return(bc)
}

bridging_coeff <- calc_bridging_coeff(g)
bridging_centrality <- bridging_coeff * betweenness(g, weights = dist_weights, normalized = FALSE)

log_msg("Calculating Vulnerability (Nodal Efficiency)...")

# Nodal Efficiency E(i) = 1/(N-1) * sum(1/d_ij)
# where d_ij is the shortest path between i and j
calc_nodal_efficiency <- function(graph) {
  # Use inverse weights for shortest path distance
  # If TOM is high, distance is small
  dists <- distances(graph, weights = 1 - E(graph)$weight)
  dists[dists == 0] <- Inf # Self-loops or same node
  
  N <- vcount(graph)
  # rowMeans of 1/dists (excluding diagonal which is Inf)
  nodal_eff <- rowSums(1/dists, na.rm = TRUE) / (N - 1)
  return(nodal_eff)
}

nodal_efficiency_val <- calc_nodal_efficiency(g)

advanced_dt <- data.table(
  taxon = names(bridging_centrality),
  bridging_centrality = bridging_centrality,
  vulnerability = nodal_efficiency_val
)

# ── 6. Merge & Annotate ───────────────────────────────────────────────────────

log_msg("Merging results and annotating...")

final_dt <- zp_dt[centrality_dt, on = "taxon"][advanced_dt, on = "taxon"]
final_dt <- mods[final_dt, on = "taxon"]

# Annotate with metadata and names
final_dt <- tax_meta[final_dt, on = "taxon"]
final_dt <- tax_map[final_dt, on = "taxon"]

# Identify likely keystone species
final_dt[, importance_rank := rank(-pagerank) + rank(-vulnerability) + rank(-betweenness)]
final_dt[, is_potential_keystone := importance_rank <= 50] 

fwrite(final_dt, file.path(STATS_OUT, "network_metrics_summary.tsv"), sep = "\t")

# ── 7. Plotting ───────────────────────────────────────────────────────────────

log_msg("Generating Z-P Plots...")

p_zp <- ggplot(final_dt, aes(x = p_score, y = z_score, color = module)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.62, linetype = "dashed", color = "red") +
  labs(title = "Network Z-P Plot", x = "Participation Coefficient (P)", y = "Within-module Degree (Z)") +
  geom_text_repel(data = final_dt[is_potential_keystone == TRUE | z_score > 3 | p_score > 0.8], 
                  aes(label = species), size = 3, max.overlaps = 15)

ggsave(file.path(STATS_OUT, "zp_plot_all.png"), p_zp, width = 10, height = 8)

p_zp_facet <- p_zp + facet_wrap(~module) + theme(legend.position = "none")
ggsave(file.path(STATS_OUT, "zp_plot_by_module.png"), p_zp_facet, width = 12, height = 10)

log_msg("Done. Results in ", STATS_OUT)
