#!/usr/bin/env Rscript
# 08_network_statistics.R — Compute network topology and centrality metrics
#
# Pipeline:
#   1. Load WGCNA consensus network and module assignments
#   2. Construct igraph object from Topological Overlap Matrix (TOM)
#   3. Calculate Within-Module Degree (z) and Participation Coefficient (p)
#   4. Calculate Centrality metrics (PageRank, Closeness, Betweenness)
#   5. Calculate Bridging Centrality and Vulnerability (Nodal Efficiency)
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

start_time <- Sys.time()

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading WGCNA results...")
net  <- readRDS("results/stage1/wgcna/consensus_wgcna.rds")
mods <- fread("results/stage1/wgcna/module_assignments.tsv")
tax_meta <- fread("results/stage1/prokaryotes_taxa_metadata.tsv")

log_msg("Extracting taxon names from raw damage summary...")
# Using subspecies as the taxon ID to match VST and mods
tax_map <- unique(fread(UPSTREAM$tax_damage, select = c("subspecies", "genus", "species")))
setnames(tax_map, "subspecies", "taxon")
tax_map <- tax_map[, .(genus = genus[1], species = species[1]), by = taxon]

vst <- readRDS("results/stage1/prokaryotes_vst.rds")
common_taxa <- intersect(colnames(vst), mods$taxon)
vst_sub <- vst[, common_taxa]

log_msg(sprintf("Processing %d taxa aligned across VST and modules", length(common_taxa)))

# ── 2. Construct Network ──────────────────────────────────────────────────────

log_msg("Calculating Adjacency and TOM...")
t1 <- Sys.time()
soft_power <- net$power
if(is.null(soft_power)) soft_power <- 12 

adj <- adjacency(vst_sub, power = soft_power, type = "signed")
tom <- TOMsimilarity(adj, TOMType = "signed")
rownames(tom) <- colnames(tom) <- common_taxa
log_msg(sprintf("  TOM calculation took %.2f seconds", as.numeric(difftime(Sys.time(), t1, units="secs"))))

# Sparse the network for efficiency
t2 <- Sys.time()
tom_sparse <- tom
tom_sparse[tom_sparse < 0.05] <- 0
g <- graph_from_adjacency_matrix(tom_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
log_msg(sprintf("  Graph construction took %.2f seconds (%d nodes, %d edges)", 
                as.numeric(difftime(Sys.time(), t2, units="secs")), vcount(g), ecount(g)))

# ── 3. Z-P Metrics ────────────────────────────────────────────────────────────

log_msg("Calculating Z and P scores...")
t3 <- Sys.time()

calc_z_p <- function(adj_mat, module_assigns) {
  taxa <- rownames(adj_mat)
  unique_mods <- unique(module_assigns$module)
  
  k_total <- rowSums(adj_mat)
  k_within <- matrix(0, nrow = length(taxa), ncol = length(unique_mods))
  colnames(k_within) <- unique_mods
  rownames(k_within) <- taxa
  
  for (m in unique_mods) {
    m_taxa <- intersect(module_assigns[module == m, taxon], taxa)
    if (length(m_taxa) > 1) {
      k_within[, m] <- rowSums(adj_mat[, m_taxa])
    } else if (length(m_taxa) == 1) {
      k_within[, m] <- adj_mat[, m_taxa]
    }
  }
  
  # Extract k_i_mi (degree within its OWN module)
  # Efficiently using match
  mod_idx <- match(taxa, module_assigns$taxon)
  node_mods <- module_assigns$module[mod_idx]
  
  ki_mi <- sapply(seq_along(taxa), function(i) {
    k_within[i, node_mods[i]]
  })
  names(ki_mi) <- taxa
  
  # Z-score within module
  z_scores <- sapply(taxa, function(x) {
    m <- node_mods[which(taxa == x)]
    m_taxa <- intersect(module_assigns[module == m, taxon], taxa)
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
log_msg(sprintf("  Z-P calculation took %.2f seconds", as.numeric(difftime(Sys.time(), t3, units="secs"))))

# ── 4. Centrality Metrics ─────────────────────────────────────────────────────

log_msg("Calculating Centrality (PageRank, Closeness, Betweenness)...")
t4 <- Sys.time()

# Distance weights for path-based metrics
dist_weights <- 1 - E(g)$weight
closeness_val <- closeness(g, weights = dist_weights, normalized = TRUE)
betweenness_val <- betweenness(g, weights = dist_weights, normalized = TRUE)
pagerank_val <- page_rank(g, weights = E(g)$weight)$vector

centrality_dt <- data.table(
  taxon = V(g)$name,
  closeness = closeness_val,
  betweenness = betweenness_val,
  pagerank = pagerank_val
)
log_msg(sprintf("  Centrality calculation took %.2f seconds", as.numeric(difftime(Sys.time(), t4, units="secs"))))

# ── 5. Bridging & Vulnerability ───────────────────────────────────────────────

log_msg("Calculating Bridging Centrality...")
t5 <- Sys.time()

calc_bridging_coeff <- function(graph) {
  deg <- degree(graph)
  bc <- sapply(V(graph), function(v) {
    neighs <- neighbors(graph, v)
    if (length(neighs) == 0) return(0)
    # Bridging coefficient measures how nodes are connected to different neighborhoods
    (1/deg[v]) / sum(1/deg[neighs])
  })
  return(bc)
}

bridging_coeff <- calc_bridging_coeff(g)
bridging_centrality <- bridging_coeff * betweenness(g, weights = dist_weights, normalized = FALSE)
log_msg(sprintf("  Bridging Centrality took %.2f seconds", as.numeric(difftime(Sys.time(), t5, units="secs"))))

log_msg("Calculating Vulnerability (Nodal Efficiency)...")
t6 <- Sys.time()

calc_nodal_efficiency <- function(graph) {
  # Efficiency is the average of inverse distances
  dists <- distances(graph, weights = 1 - E(graph)$weight)
  dists[dists == 0] <- Inf 
  N <- vcount(graph)
  nodal_eff <- rowSums(1/dists, na.rm = TRUE) / (N - 1)
  return(nodal_eff)
}

nodal_efficiency_val <- calc_nodal_efficiency(g)
advanced_dt <- data.table(
  taxon = V(g)$name,
  bridging_centrality = bridging_centrality,
  vulnerability = nodal_efficiency_val
)
log_msg(sprintf("  Nodal Efficiency took %.2f seconds", as.numeric(difftime(Sys.time(), t6, units="secs"))))

# ── 6. Merge & Annotate ───────────────────────────────────────────────────────

log_msg("Merging results and annotating...")
t7 <- Sys.time()

# Use merge() for safety and to avoid duplicate column names
final_dt <- merge(zp_dt, centrality_dt, by = "taxon", all = TRUE)
final_dt <- merge(final_dt, advanced_dt, by = "taxon", all = TRUE)
final_dt <- merge(final_dt, mods, by = "taxon", all.x = TRUE)
final_dt <- merge(final_dt, tax_meta, by = "taxon", all.x = TRUE)
final_dt <- merge(final_dt, tax_map, by = "taxon", all.x = TRUE)

# Handle potential NAs in metrics for ranking
final_dt[is.na(pagerank), pagerank := 0]
final_dt[is.na(vulnerability), vulnerability := 0]
final_dt[is.na(betweenness), betweenness := 0]
final_dt[is.na(closeness), closeness := 0]

# 1. Identify Keystones: Top 50 by sum of ranks (PageRank + Vulnerability + Betweenness)
final_dt[, rank_sum := rank(-pagerank) + rank(-vulnerability) + rank(-betweenness)]
final_dt[, is_potential_keystone := rank(rank_sum) <= 50]

# 2. Identify Hidden Gems: High influence (PageRank/Closeness) but low hub-ness (Z-score)
# Let's say top 5% PageRank but Z-score < 1.0
pr_threshold <- quantile(final_dt$pagerank, 0.95)
final_dt[, is_hidden_gem := pagerank >= pr_threshold & z_score < 1.0]

# 3. Categorize importance for easy filtering
final_dt[, importance_type := fcase(
  is_potential_keystone, "Keystone",
  role %in% c("Module Hub", "Network Hub"), "Hub",
  is_hidden_gem, "Hidden Gem",
  role == "Connector", "Connector",
  default = "Peripheral"
)]

fwrite(final_dt, file.path(STATS_OUT, "network_metrics_summary.tsv"), sep = "\t")
log_msg(sprintf("  Merging and saving took %.2f seconds", as.numeric(difftime(Sys.time(), t7, units="secs"))))

# ── 7. Plotting ───────────────────────────────────────────────────────────────

log_msg("Generating Z-P Plots...")

p_zp <- ggplot(final_dt, aes(x = p_score, y = z_score, color = module)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.62, linetype = "dashed", color = "red") +
  labs(title = "Network Z-P Plot", x = "Participation Coefficient (P)", y = "Within-module Degree (Z)") +
  geom_text_repel(data = final_dt[(is_potential_keystone == TRUE | z_score > 3 | p_score > 0.8) & !is.na(species)], 
                  aes(label = species), size = 3, max.overlaps = 20)

ggsave(file.path(STATS_OUT, "zp_plot_all.png"), p_zp, width = 10, height = 8)

p_zp_facet <- p_zp + facet_wrap(~module) + theme(legend.position = "none")
ggsave(file.path(STATS_OUT, "zp_plot_by_module.png"), p_zp_facet, width = 12, height = 10)

total_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
log_msg(sprintf("Done. Total runtime: %.2f minutes. Results in %s", total_time, STATS_OUT))
