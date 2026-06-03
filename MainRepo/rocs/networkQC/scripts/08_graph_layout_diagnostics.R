#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
  library(igraph)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()

OUT_FIG <- here("networkQC", "results", "figures", "layout_diagnostics")
OUT_TABLE <- here("networkQC", "results", "tables")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
train <- do.call(rbind, lapply(PARAMS$stage1_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
}))

settings <- data.table(
  setting = c("baseline", "opt5", "leiden"),
  power = c(20, 12, 12),
  assign_path = c(
    here("networkQC", "results", "full_eval", "baseline", "module_assignments.tsv"),
    here("networkQC", "results", "full_eval", "opt5", "module_assignments.tsv"),
    here("networkQC", "results", "tables", "leiden_module_assignments.tsv")
  )
)

module_palette <- function(mods) {
  known <- c(
    grey = "grey80", turquoise = "turquoise3", blue = "royalblue3",
    brown = "sienna4", yellow = "gold2", green = "forestgreen"
  )
  lev <- sort(unique(mods))
  out <- setNames(rainbow(length(lev)), lev)
  hit <- intersect(names(known), lev)
  out[hit] <- known[hit]
  out
}

build_graph <- function(power, assign_path, edge_q = 0.995) {
  mods <- fread(assign_path)
  modmap <- setNames(mods$module, mods$taxon)
  adj <- adjacency(train, power = power, type = "signed")
  tom <- TOMsimilarity(adj, TOMType = "signed")
  diag(tom) <- 0
  thr <- quantile(tom[upper.tri(tom)], edge_q, na.rm = TRUE)
  idx <- which(tom >= thr, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  edges <- data.table(
    from = colnames(train)[idx[, 1]],
    to = colnames(train)[idx[, 2]],
    weight = tom[idx]
  )
  edges[, from_module := modmap[from]]
  edges[, to_module := modmap[to]]
  edges[, within_module := from_module == to_module]
  g <- graph_from_data_frame(edges[, .(from, to, weight)], directed = FALSE,
                             vertices = data.frame(name = colnames(train)))
  V(g)$module <- modmap[V(g)$name]
  E(g)$weight <- edges$weight
  list(g = g, edges = edges, threshold = thr)
}

plot_graph <- function(g, file, title, layout_kind = "fr") {
  pal <- module_palette(V(g)$module)
  vcol <- pal[V(g)$module]
  set.seed(PARAMS$seed)
  lay <- switch(
    layout_kind,
    kk = layout_with_kk(g, weights = E(g)$weight),
    circle = layout_in_circle(g),
    layout_with_fr(g, weights = E(g)$weight, niter = 2000)
  )
  png(file, width = 1600, height = 1200, res = 150)
  plot(g, layout = lay, vertex.size = 3, vertex.label = NA,
       vertex.color = vcol, vertex.frame.color = NA,
       edge.width = pmax(0.2, E(g)$weight / max(E(g)$weight, na.rm = TRUE) * 1.2),
       edge.color = rgb(0, 0, 0, 0.12),
       main = title)
  legend("topleft", legend = names(pal), col = pal, pch = 16, pt.cex = 1.2,
         cex = 0.8, bty = "n")
  dev.off()
}

plot_module_graph <- function(edges, file, title) {
  mg <- edges[, .(weight = sum(weight), n_edges = .N), by = .(from_module, to_module)]
  mg[, a := pmin(from_module, to_module)]
  mg[, b := pmax(from_module, to_module)]
  mg <- mg[, .(weight = sum(weight), n_edges = sum(n_edges)), by = .(a, b)]
  g <- graph_from_data_frame(mg[, .(from = a, to = b, weight)], directed = FALSE)
  pal <- module_palette(V(g)$name)
  png(file, width = 1200, height = 900, res = 150)
  set.seed(PARAMS$seed)
  plot(g, layout = layout_with_fr(g, weights = E(g)$weight, niter = 2000),
       vertex.size = 22, vertex.label.cex = 0.85,
       vertex.color = pal[V(g)$name], vertex.frame.color = "white",
       edge.width = pmax(1, E(g)$weight / max(E(g)$weight) * 8),
       edge.color = rgb(0, 0, 0, 0.22), main = title)
  dev.off()
}

diagnostics <- rbindlist(lapply(seq_len(nrow(settings)), function(i) {
  s <- settings[i]
  bg <- build_graph(s$power, s$assign_path)
  g <- bg$g
  edges <- bg$edges
  deg <- degree(g)
  g_noniso <- induced_subgraph(g, vids = V(g)[deg > 0])
  comp <- components(g_noniso)
  largest_id <- which.max(comp$csize)
  g_largest <- induced_subgraph(g_noniso, vids = V(g_noniso)[comp$membership == largest_id])
  g_nongrey <- induced_subgraph(g_noniso, vids = V(g_noniso)[module != "grey"])

  plot_graph(g_noniso, file.path(OUT_FIG, paste0(s$setting, "_nonisolates_fr.png")),
             paste0(s$setting, ": non-isolates, FR layout"), "fr")
  plot_graph(g_noniso, file.path(OUT_FIG, paste0(s$setting, "_nonisolates_kk.png")),
             paste0(s$setting, ": non-isolates, KK layout"), "kk")
  plot_graph(g_largest, file.path(OUT_FIG, paste0(s$setting, "_largest_component_fr.png")),
             paste0(s$setting, ": largest connected component"), "fr")
  if (vcount(g_nongrey) > 2 && ecount(g_nongrey) > 0) {
    plot_graph(g_nongrey, file.path(OUT_FIG, paste0(s$setting, "_nongrey_nonisolates_fr.png")),
               paste0(s$setting, ": non-grey non-isolates"), "fr")
  }
  plot_module_graph(edges, file.path(OUT_FIG, paste0(s$setting, "_module_meta_graph.png")),
                    paste0(s$setting, ": module-level edge graph"))

  data.table(
    setting = s$setting,
    power = s$power,
    edge_threshold = bg$threshold,
    nodes_total = vcount(g),
    edges = ecount(g),
    nodes_nonisolated = vcount(g_noniso),
    isolates = sum(deg == 0),
    nonisolated_pct = vcount(g_noniso) / vcount(g) * 100,
    components_nonisolated = components(g_noniso)$no,
    largest_component_nodes = vcount(g_largest),
    within_module_edge_pct = mean(edges$within_module, na.rm = TRUE) * 100
  )
}), fill = TRUE)

fwrite(diagnostics, file.path(OUT_TABLE, "graph_layout_diagnostics.tsv"), sep = "\t")
message("Graph layout diagnostics written to ", OUT_FIG)

