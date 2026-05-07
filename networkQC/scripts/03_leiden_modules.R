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

OUT_TABLE <- here("networkQC", "results", "tables")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

if (!"cluster_leiden" %in% getNamespaceExports("igraph")) {
  stop("igraph::cluster_leiden not available in this environment.")
}

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
wgcna_mods <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores
train <- do.call(rbind, lapply(PARAMS$stage1_cores, function(c) expr_by_core[[c]]))

reco_file <- file.path(OUT_TABLE, "qc_parameter_sweep_recommended.tsv")
power <- 20
if (file.exists(reco_file)) {
  reco <- fread(reco_file)
  if (nrow(reco) > 0 && is.finite(reco$power[1])) power <- as.integer(reco$power[1])
}
log_msg(sprintf("Using power=%d for Leiden network backbone", power))

tom <- TOMsimilarityFromExpr(train, networkType = "signed", TOMType = "signed", power = power, verbose = 0)
diag(tom) <- 0
thr <- 0.05
idx <- which(tom >= thr, arr.ind = TRUE)
idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
edges <- data.table(
  from = colnames(train)[idx[, 1]],
  to = colnames(train)[idx[, 2]],
  weight = tom[idx]
)
g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = colnames(train)))

resolutions <- c(0.05, 0.10, 0.20, 0.40)
leiden_summary <- rbindlist(lapply(resolutions, function(r) {
  cl <- cluster_leiden(g, objective_function = "modularity", resolution = r, weights = E(g)$weight)
  mem <- membership(cl)
  tab <- table(mem)
  data.table(
    resolution = r,
    n_modules = length(tab),
    module_size_median = median(as.numeric(tab)),
    largest_module = max(as.numeric(tab))
  )
}))
fwrite(leiden_summary, file.path(OUT_TABLE, "leiden_resolution_summary.tsv"), sep = "\t")

best_res <- leiden_summary[order(abs(n_modules - wgcna_mods[module != "grey", uniqueN(module)]), -module_size_median)][1, resolution]
cl <- cluster_leiden(g, objective_function = "modularity", resolution = best_res, weights = E(g)$weight)
mem <- membership(cl)
leid_mods <- data.table(
  taxon = names(mem),
  module = paste0("L", mem)
)
fwrite(leid_mods, file.path(OUT_TABLE, "leiden_module_assignments.tsv"), sep = "\t")

meta_s <- meta[, .(label, core, age_kyr = y_bp / 1000)]
all_samples <- do.call(rbind, expr_by_core)
me_leid <- moduleEigengenes(all_samples, colors = setNames(leid_mods$module, leid_mods$taxon)[colnames(all_samples)])$eigengenes
me_leid <- as.data.table(me_leid, keep.rownames = "sample")
fwrite(me_leid, file.path(OUT_TABLE, "leiden_module_eigengenes.tsv"), sep = "\t")

summary_dt <- data.table(
  metric = c("selected_power", "tom_threshold", "selected_resolution", "leiden_modules"),
  value = c(power, thr, best_res, uniqueN(leid_mods$module))
)
fwrite(summary_dt, file.path(OUT_TABLE, "leiden_run_summary.tsv"), sep = "\t")
log_msg("Leiden module construction complete")

