#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Resolution tuning requires igraph (cluster_leiden backend).")
}
if (!("cluster_leiden" %in% getNamespaceExports("igraph"))) {
  stop("Resolution tuning requires igraph::cluster_leiden, but it is unavailable in this R environment.")
}

sweep_dir <- file.path(DIRS$results, "resolution_sweep")
ensure_dirs(c(sweep_dir))

meta_path <- file.path(DIRS$main, "sample_metadata_main.tsv")
clr_path <- file.path(DIRS$main, "clr_matrix_train_centered.rds")
if (!file.exists(meta_path) || !file.exists(clr_path)) {
  stop(
    "Missing main-stage inputs for sweep. Expected files:\n",
    " - ", meta_path, "\n",
    " - ", clr_path, "\n",
    "Run Leiden 01_data_prep first for this results root."
  )
}

prev_meta_path <- file.path(DIRS$main, "02_wgcna_main_run_metadata.tsv")
prev_meta <- if (file.exists(prev_meta_path)) fread(prev_meta_path) else data.table()
prev_resolution_requested <- if (nrow(prev_meta) > 0 && "leiden_resolution" %in% names(prev_meta)) prev_meta$leiden_resolution[1] else NA_character_
prev_resolution_used <- if (nrow(prev_meta) > 0 && "leiden_resolution_used" %in% names(prev_meta)) prev_meta$leiden_resolution_used[1] else NA_character_

meta_main <- fread(meta_path)
clr <- readRDS(clr_path)

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  ids <- meta_main[core == core_id, label]
  ids <- intersect(ids, rownames(clr))
  clr[ids, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

powers <- c(1:10, seq(12, 20, by = 2))
sft_list <- lapply(PARAMS$training_cores, function(core_id) {
  sft <- pickSoftThreshold(expr_by_core[[core_id]], powerVector = powers, networkType = "signed", verbose = 0)
  dt <- as.data.table(sft$fitIndices)
  dt[, core := core_id]
  dt
})
sft_dt <- rbindlist(sft_list)
sft_dt[, signedR2 := ifelse(slope < 0, SFT.R.sq, 0)]

passing <- sft_dt[signedR2 >= PARAMS$soft_power_target_r2, .(min_pass = min(Power)), by = core]
if (nrow(passing) == length(PARAMS$training_cores)) {
  soft_power <- max(passing$min_pass)
} else {
  best <- sft_dt[slope < 0 & mean.k. > 5, .(best_r2 = max(signedR2)), by = Power]
  soft_power <- if (nrow(best) > 0) best[which.max(best_r2), Power] else 12L
}

multiExpr <- lapply(PARAMS$training_cores, function(core_id) list(data = expr_by_core[[core_id]]))
names(multiExpr) <- PARAMS$training_cores
train_expr <- do.call(rbind, lapply(PARAMS$training_cores, function(c) expr_by_core[[c]]))

cons <- consensusDissTOMandTree(multiExpr = multiExpr, softPower = soft_power)
consensus_matrix <- cons$consensusTOM
if (is.null(consensus_matrix)) {
  stop("WGCNA consensusDissTOMandTree did not return consensusTOM; cannot run Leiden sweep.")
}
consensus_matrix[!is.finite(consensus_matrix)] <- 0
consensus_matrix <- (consensus_matrix + t(consensus_matrix)) / 2
diag(consensus_matrix) <- 0

matrix_mode <- PARAMS$leiden_consensus_matrix_interpretation
if (identical(matrix_mode, "auto")) {
  med_val <- stats::median(consensus_matrix[upper.tri(consensus_matrix)], na.rm = TRUE)
  matrix_mode <- if (is.finite(med_val) && med_val > 0.5) "dissimilarity" else "similarity"
}
if (identical(matrix_mode, "dissimilarity")) {
  consensus_sim <- 1 - consensus_matrix
} else if (identical(matrix_mode, "similarity")) {
  consensus_sim <- consensus_matrix
} else {
  stop("Unsupported consensus matrix interpretation: ", matrix_mode)
}
consensus_sim[!is.finite(consensus_sim)] <- 0
consensus_sim[consensus_sim < 0] <- 0
consensus_sim[consensus_sim > 1] <- 1
diag(consensus_sim) <- 0

k_neighbors <- min(as.integer(PARAMS$leiden_graph_k_neighbors), ncol(consensus_sim) - 1L)
sparse_adj <- matrix(0, nrow = nrow(consensus_sim), ncol = ncol(consensus_sim), dimnames = dimnames(consensus_sim))
for (i in seq_len(nrow(consensus_sim))) {
  ord <- order(consensus_sim[i, ], decreasing = TRUE)
  ord <- ord[ord != i]
  keep <- ord[seq_len(k_neighbors)]
  sparse_adj[i, keep] <- consensus_sim[i, keep]
}
graph_sym <- PARAMS$leiden_graph_symmetrization
if (identical(graph_sym, "union")) {
  sparse_adj <- pmax(sparse_adj, t(sparse_adj))
} else if (identical(graph_sym, "mutual")) {
  sparse_adj <- ifelse((sparse_adj > 0) & (t(sparse_adj) > 0), consensus_sim, 0)
} else {
  stop("Unsupported Leiden graph symmetrization: ", graph_sym)
}
diag(sparse_adj) <- 0

g <- igraph::graph_from_adjacency_matrix(sparse_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
leiden_formals <- names(formals(igraph::cluster_leiden))

resolution_grid <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0)
resolution_env <- env_csv("WGCNA_HMM_LEIDEN_SWEEP_GRID", default = character())
if (length(resolution_env) > 0) {
  resolution_grid <- suppressWarnings(as.numeric(resolution_env))
}
resolution_grid <- sort(unique(resolution_grid[is.finite(resolution_grid) & resolution_grid > 0]))
if (length(resolution_grid) == 0) stop("No valid positive resolutions to test.")

run_leiden_summary <- function(resolution_value) {
  leiden_args <- list(graph = g, weights = igraph::E(g)$weight)
  if ("objective_function" %in% leiden_formals) leiden_args$objective_function <- "modularity"
  if ("resolution" %in% leiden_formals) {
    leiden_args$resolution <- resolution_value
  } else if ("resolution_parameter" %in% leiden_formals) {
    leiden_args$resolution_parameter <- resolution_value
  }
  if ("n_iterations" %in% leiden_formals) leiden_args$n_iterations <- 100

  cl <- do.call(igraph::cluster_leiden, leiden_args)
  memb <- igraph::membership(cl)
  memb_dt <- as.data.table(table(memb))
  setnames(memb_dt, c("community", "n_taxa"))
  memb_dt[, community := as.character(community)]
  n_communities_raw <- nrow(memb_dt)
  kept_dt <- memb_dt[n_taxa >= PARAMS$wgcna_min_module_size]
  n_modules_non_grey <- nrow(kept_dt)
  non_grey_taxa <- kept_dt[, sum(n_taxa)]
  if (!is.finite(non_grey_taxa)) non_grey_taxa <- 0L
  grey_size <- length(memb) - non_grey_taxa
  grey_fraction <- grey_size / length(memb)

  non_grey_sizes <- kept_dt$n_taxa
  min_non_grey <- if (length(non_grey_sizes) > 0) min(non_grey_sizes) else NA_real_
  median_non_grey <- if (length(non_grey_sizes) > 0) stats::median(non_grey_sizes) else NA_real_
  max_non_grey <- if (length(non_grey_sizes) > 0) max(non_grey_sizes) else NA_real_
  largest_non_grey_fraction <- if (length(non_grey_sizes) > 0) max(non_grey_sizes) / length(memb) else NA_real_

  modularity_score <- suppressWarnings(igraph::modularity(g, membership = memb, weights = igraph::E(g)$weight))

  data.table(
    resolution = resolution_value,
    n_taxa = length(memb),
    n_communities_raw = n_communities_raw,
    n_modules_non_grey = n_modules_non_grey,
    n_modules_surviving_min_size = n_modules_non_grey,
    non_grey_taxa = non_grey_taxa,
    grey_size = grey_size,
    grey_fraction = grey_fraction,
    min_non_grey_size = min_non_grey,
    median_non_grey_size = median_non_grey,
    max_non_grey_size = max_non_grey,
    largest_non_grey_fraction = largest_non_grey_fraction,
    modularity = modularity_score,
    eigengenes_proxy_n = n_modules_non_grey + 1L,
    hmm_runnable_proxy = n_modules_non_grey >= 3
  )
}

sweep_dt <- rbindlist(lapply(resolution_grid, run_leiden_summary), use.names = TRUE)
setorder(sweep_dt, resolution)

sweep_dt[, feasible :=
  hmm_runnable_proxy == TRUE &
  n_modules_non_grey >= 3 &
  n_modules_non_grey <= 30 &
  grey_fraction <= 0.70 &
  (is.na(largest_non_grey_fraction) | largest_non_grey_fraction <= 0.60)
]

feasible_dt <- sweep_dt[feasible == TRUE]
selection_rule <- c(
  "1) Candidate must pass guardrails: non-grey modules in [3,30], grey fraction <= 0.70, largest non-grey module <= 60% of taxa, HMM proxy runnable.",
  "2) Among feasible candidates, retain those within 0.01 modularity of the best feasible modularity.",
  "3) From that near-optimal set, choose the fewest non-grey modules (interpretability).",
  "4) Tie-break by higher modularity, then lower resolution."
)

if (nrow(feasible_dt) > 0) {
  max_mod <- max(feasible_dt$modularity, na.rm = TRUE)
  near_opt <- feasible_dt[modularity >= (max_mod - 0.01)]
  setorder(near_opt, n_modules_non_grey, -modularity, resolution)
  selected <- near_opt[1]
  selection_status <- "selected_from_feasible"
} else {
  fallback <- sweep_dt[hmm_runnable_proxy == TRUE]
  if (nrow(fallback) == 0) fallback <- copy(sweep_dt)
  setorder(fallback, -modularity, n_modules_non_grey, resolution)
  selected <- fallback[1]
  selection_status <- "fallback_best_modularity_no_feasible"
}

sweep_dt[, selected := resolution == selected$resolution]

fwrite(sweep_dt, file.path(sweep_dir, "leiden_resolution_sweep.tsv"), sep = "\t")
fwrite(selected, file.path(sweep_dir, "leiden_resolution_selected.tsv"), sep = "\t")

summary_lines <- c(
  "# Leiden resolution sweep selection summary",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Provenance",
  sprintf("- Results root: `%s`", DIRS$results),
  sprintf("- Script: `%s`", "code/wgcna_hmm/10_tune_leiden_resolution.R"),
  sprintf("- Training cores: %s", paste(PARAMS$training_cores, collapse = ", ")),
  sprintf("- Min module size threshold: %d", PARAMS$wgcna_min_module_size),
  sprintf("- Consensus graph: %s kNN sparse graph from interpreted consensus matrix (k=%d)", graph_sym, k_neighbors),
  sprintf("- Consensus matrix interpretation: %s", matrix_mode),
  sprintf("- Previous run metadata requested resolution: %s", ifelse(is.na(prev_resolution_requested), "NA", prev_resolution_requested)),
  sprintf("- Previous run metadata used resolution in practice: %s", ifelse(is.na(prev_resolution_used), "NA", prev_resolution_used)),
  "",
  "## Resolutions tested",
  paste0("- ", paste(format(resolution_grid, trim = TRUE, scientific = FALSE), collapse = ", ")),
  "",
  "## Selection rule (conservative and explicit)",
  paste0("- ", selection_rule),
  "",
  "## Selected resolution",
  sprintf("- Selected resolution: **%.3f**", selected$resolution),
  sprintf("- Selection status: `%s`", selection_status),
  sprintf("- Modularity: %.5f", selected$modularity),
  sprintf("- Non-grey modules: %d", selected$n_modules_non_grey),
  sprintf("- Grey fraction: %.4f", selected$grey_fraction),
  sprintf("- Non-grey size summary (min/median/max): %s / %s / %s",
          ifelse(is.na(selected$min_non_grey_size), "NA", as.character(selected$min_non_grey_size)),
          ifelse(is.na(selected$median_non_grey_size), "NA", as.character(selected$median_non_grey_size)),
          ifelse(is.na(selected$max_non_grey_size), "NA", as.character(selected$max_non_grey_size))),
  "",
  "## Notes",
  "- This is a sensitivity/alternative partition workflow; baseline results/wgcna_hmm remain unchanged.",
  "- Full sweep table is in `leiden_resolution_sweep.tsv`."
)

writeLines(summary_lines, con = file.path(sweep_dir, "leiden_resolution_selection_summary.md"))

write_run_metadata(
  file.path(sweep_dir, "10_tune_leiden_resolution_run_metadata.tsv"),
  "10_tune_leiden_resolution.R",
  extra = list(
    results_root = DIRS$results,
    soft_power = soft_power,
    leiden_graph_k_neighbors = k_neighbors,
    leiden_graph_symmetrization = graph_sym,
    leiden_consensus_matrix_interpretation = matrix_mode,
    resolutions_tested = paste(resolution_grid, collapse = ";"),
    selected_resolution = selected$resolution,
    selection_status = selection_status,
    previous_requested_resolution = prev_resolution_requested,
    previous_used_resolution = prev_resolution_used
  )
)
write_session_info(file.path(sweep_dir, "10_tune_leiden_resolution_sessionInfo.txt"))

log_msg("10_tune_leiden_resolution complete. Selected resolution=", selected$resolution)
