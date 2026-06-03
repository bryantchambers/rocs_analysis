#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(tidyverse)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_mode <- PARAMS$wgcna_run_mode
if (length(args) > 0) {
  mode_arg <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
  if (!is.na(mode_arg) && nzchar(mode_arg)) run_mode <- mode_arg
}
if (!run_mode %in% c("build", "final")) {
  stop("Invalid WGCNA mode: ", run_mode, ". Use --mode=build or --mode=final.")
}
preservation_perms <- if (run_mode == "final") {
  PARAMS$wgcna_preservation_permutations_final
} else {
  PARAMS$wgcna_preservation_permutations_build
}
log_msg(sprintf("WGCNA mode: %s; preservation permutations: %d", run_mode, preservation_perms))

clr <- readRDS(file.path(DIRS$main, "clr_matrix_train_centered.rds"))
meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  ids <- meta_main[core == core_id, label]
  ids <- intersect(ids, rownames(clr))
  clr[ids, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

training_manifest_path <- file.path(DIRS$main, "wgcna_training_samples.tsv")
if (!file.exists(training_manifest_path)) {
  stop("Missing WGCNA training manifest: ", training_manifest_path, ". Run 01_data_prep.R first.")
}
wgcna_training_manifest <- fread(training_manifest_path)
required_manifest_cols <- c("label", "core")
missing_manifest_cols <- setdiff(required_manifest_cols, names(wgcna_training_manifest))
if (length(missing_manifest_cols) > 0) {
  stop("WGCNA training manifest missing columns: ", paste(missing_manifest_cols, collapse = ", "))
}

fit_expr_by_core <- lapply(PARAMS$training_cores, function(core_id) {
  ids <- wgcna_training_manifest[core == core_id, label]
  ids <- intersect(ids, rownames(clr))
  if (length(ids) == 0) stop("No selected WGCNA training samples for core: ", core_id)
  clr[ids, , drop = FALSE]
})
names(fit_expr_by_core) <- PARAMS$training_cores
fit_core_n <- vapply(fit_expr_by_core, nrow, integer(1))
if (identical(PARAMS$wgcna_input_strategy, "balanced") && length(unique(fit_core_n)) != 1L) {
  stop("Balanced WGCNA manifest does not have equal selected counts across training cores.")
}
log_msg(
  "WGCNA input strategy: ", PARAMS$wgcna_input_strategy,
  "; profile: ", PARAMS$wgcna_profile,
  "; selected training counts: ", paste(names(fit_core_n), fit_core_n, sep = "=", collapse = ", "),
  "; HMM input balancing: ", PARAMS$hmm_input_balancing_status
)

powers <- c(1:10, seq(12, 20, by = 2))
sft_list <- lapply(PARAMS$training_cores, function(core_id) {
  sft <- pickSoftThreshold(fit_expr_by_core[[core_id]], powerVector = powers, networkType = "signed", verbose = 0)
  dt <- as.data.table(sft$fitIndices)
  dt[, core := core_id]
  dt
})
sft_dt <- rbindlist(sft_list)
sft_dt[, signedR2 := ifelse(slope < 0, SFT.R.sq, 0)]

if (!is.null(PARAMS$wgcna_soft_power) && is.finite(PARAMS$wgcna_soft_power)) {
  soft_power <- as.integer(PARAMS$wgcna_soft_power)
  soft_power_selection <- "fixed_config"
  log_msg(sprintf("Using fixed WGCNA soft power from config: %d", soft_power))
} else {
  passing <- sft_dt[signedR2 >= PARAMS$soft_power_target_r2, .(min_pass = min(Power)), by = core]
  if (nrow(passing) == length(PARAMS$training_cores)) {
    soft_power <- max(passing$min_pass)
    soft_power_selection <- "all_training_cores_pass_target"
  } else {
    best <- sft_dt[slope < 0 & mean.k. > 5, .(best_r2 = max(signedR2)), by = Power]
    soft_power <- if (nrow(best) > 0) best[which.max(best_r2), Power] else 12L
    soft_power_selection <- "fallback_best_signed_r2_mean_k_gt_5"
  }
  log_msg(sprintf("Auto-selected WGCNA soft power: %d (%s)", soft_power, soft_power_selection))
}

multiExpr <- lapply(PARAMS$training_cores, function(core_id) list(data = fit_expr_by_core[[core_id]]))
names(multiExpr) <- PARAMS$training_cores

train_expr <- do.call(rbind, lapply(PARAMS$training_cores, function(c) fit_expr_by_core[[c]]))

if (identical(PARAMS$module_method, "wgcna")) {
  net <- blockwiseConsensusModules(
    multiExpr = multiExpr,
    power = soft_power,
    networkType = "signed",
    corType = "pearson",
    maxBlockSize = 5000,
    minModuleSize = PARAMS$wgcna_min_module_size,
    deepSplit = PARAMS$wgcna_deep_split,
    mergeCutHeight = PARAMS$wgcna_merge_cut_height,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    verbose = 2
  )
  module_colors <- net$colors
  leiden_backend <- NA_character_
} else if (identical(PARAMS$module_method, "leiden")) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Leiden module method requested but igraph is not installed in this R environment.")
  }
  if (!("cluster_leiden" %in% getNamespaceExports("igraph"))) {
    stop("Leiden module method requested but igraph::cluster_leiden is unavailable in this R environment.")
  }

  # Build a WGCNA consensus network from the same training-core expressions,
  # then run Leiden only for module partitioning of that consensus graph.
  cons <- consensusDissTOMandTree(multiExpr = multiExpr, softPower = soft_power)
  consensus_matrix <- cons$consensusTOM
  if (is.null(consensus_matrix)) {
    stop("WGCNA consensusDissTOMandTree did not return consensusTOM; cannot run Leiden on consensus network.")
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

  # Build a sparse graph from consensus TOM to avoid complete-graph collapse in modularity optimization.
  k_neighbors <- min(as.integer(PARAMS$leiden_graph_k_neighbors), ncol(consensus_sim) - 1L)
  sparse_adj <- matrix(0, nrow = nrow(consensus_sim), ncol = ncol(consensus_sim),
                       dimnames = dimnames(consensus_sim))
  for (i in seq_len(nrow(consensus_sim))) {
    ord <- order(consensus_sim[i, ], decreasing = TRUE)
    ord <- ord[ord != i]
    keep <- ord[seq_len(k_neighbors)]
    sparse_adj[i, keep] <- consensus_sim[i, keep]
  }
  # Symmetrize either by union-kNN (max weight) or true mutual-kNN (intersection of retained neighbors).
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
  run_leiden <- function(resolution_value) {
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
    memb_tab <- as.data.table(table(memb))
    setnames(memb_tab, c("community", "n_taxa"))
    memb_tab[, community := as.character(community)]
    kept <- memb_tab[n_taxa >= PARAMS$wgcna_min_module_size, community]
    list(membership = memb, membership_table = memb_tab, kept_communities = kept)
  }

  requested_resolution <- PARAMS$leiden_resolution
  chosen_resolution <- requested_resolution
  auto_adjusted <- FALSE

  if (!isTRUE(PARAMS$leiden_resolution_strict)) {
    resolution_candidates <- unique(c(
      requested_resolution,
      requested_resolution * 1.5,
      requested_resolution * 2,
      requested_resolution * 3
    ))
    leiden_trials <- lapply(resolution_candidates, run_leiden)
    n_non_grey_by_trial <- vapply(leiden_trials, function(x) length(x$kept_communities), integer(1))

    chosen_idx <- which(n_non_grey_by_trial >= 3)[1]
    if (is.na(chosen_idx)) chosen_idx <- which.max(n_non_grey_by_trial)

    chosen_resolution <- resolution_candidates[chosen_idx]
    chosen <- leiden_trials[[chosen_idx]]
    auto_adjusted <- !identical(chosen_resolution, requested_resolution)
  } else {
    chosen <- run_leiden(requested_resolution)
    n_non_grey_by_trial <- integer(0)
    resolution_candidates <- numeric(0)
  }

  memb <- chosen$membership
  memb_tab <- chosen$membership_table
  kept_communities <- chosen$kept_communities
  small_communities <- setdiff(memb_tab$community, kept_communities)

  if (auto_adjusted) {
    log_msg(
      "Leiden resolution adjusted from ", requested_resolution,
      " to ", chosen_resolution,
      " to satisfy downstream minimum module count (>=3 non-grey modules)."
    )
  }

  module_labels <- rep("grey", length(memb))
  names(module_labels) <- names(memb)
  if (length(kept_communities) > 0) {
    kept_idx <- match(as.character(memb), kept_communities)
    module_labels[!is.na(kept_idx)] <- paste0("leiden", kept_idx[!is.na(kept_idx)])
  }

  merge_threshold <- PARAMS$leiden_merge_eigengene_cor
  merge_applied <- FALSE
  merge_n_steps <- 0L
  if (is.finite(merge_threshold)) {
    if (merge_threshold <= 0 || merge_threshold > 1) {
      stop("WGCNA_HMM_LEIDEN_MERGE_EIGENGENE_COR must be in (0, 1].")
    }
    repeat {
      non_grey_modules <- sort(unique(module_labels[module_labels != "grey"]))
      if (length(non_grey_modules) < 2) break

      MEs_tmp <- moduleEigengenes(train_expr, module_labels)$eigengenes
      if (ncol(MEs_tmp) < 2) break

      me_cor <- suppressWarnings(cor(MEs_tmp, use = "pairwise.complete.obs"))
      diag(me_cor) <- NA_real_
      max_cor <- suppressWarnings(max(me_cor, na.rm = TRUE))
      if (!is.finite(max_cor) || max_cor < merge_threshold) break

      idx <- which(me_cor == max_cor, arr.ind = TRUE)[1, ]
      me_i <- colnames(me_cor)[idx[1]]
      me_j <- colnames(me_cor)[idx[2]]
      mod_i <- sub("^ME", "", me_i)
      mod_j <- sub("^ME", "", me_j)
      if (!identical(mod_i, mod_j)) {
        module_labels[module_labels == mod_j] <- mod_i
        merge_applied <- TRUE
        merge_n_steps <- merge_n_steps + 1L
      }
      if (merge_n_steps >= as.integer(PARAMS$leiden_merge_max_iter)) {
        warning("Reached WGCNA_HMM_LEIDEN_MERGE_MAX_ITER; stopping module merge loop.")
        break
      }
    }

    post_merge_tab <- as.data.table(table(module_labels))
    setnames(post_merge_tab, c("module", "n_taxa"))
    kept_after_merge <- post_merge_tab[module != "grey" & n_taxa >= PARAMS$wgcna_min_module_size, module]
    module_labels[!(module_labels %in% c("grey", kept_after_merge))] <- "grey"
  }

  module_colors <- module_labels[colnames(train_expr)]
  names(module_colors) <- colnames(train_expr)

  net <- list(
    method = "leiden",
    backend = "igraph::cluster_leiden",
    soft_power = soft_power,
    resolution_parameter_requested = requested_resolution,
    resolution_parameter_used = chosen_resolution,
    resolution_strict = PARAMS$leiden_resolution_strict,
    resolution_auto_adjusted = auto_adjusted,
    min_module_size = PARAMS$wgcna_min_module_size,
    network_type = "wgcna_consensus_tom",
    consensus_network_constructor = "WGCNA::consensusDissTOMandTree",
    leiden_graph_representation = sprintf("%s kNN sparse graph from interpreted consensus matrix (k=%d)", graph_sym, k_neighbors),
    leiden_consensus_matrix_interpretation = matrix_mode,
    leiden_graph_k_neighbors = k_neighbors,
    leiden_graph_symmetrization = graph_sym,
    leiden_merge_eigengene_cor = if (is.finite(merge_threshold)) merge_threshold else NA_real_,
    leiden_merge_applied = merge_applied,
    leiden_merge_n_steps = merge_n_steps,
    consensus_tree = cons$consTree,
    communities_raw = memb,
    module_labels = module_colors,
    n_communities_raw = uniqueN(memb),
    n_modules_non_grey = uniqueN(module_colors[module_colors != "grey"]),
    n_modules_non_grey_by_resolution_trial = as.list(stats::setNames(n_non_grey_by_trial, as.character(resolution_candidates))),
    graph_weight_summary = list(
      n_taxa = ncol(train_expr),
      edge_weight_min = min(sparse_adj[upper.tri(sparse_adj)][sparse_adj[upper.tri(sparse_adj)] > 0], na.rm = TRUE),
      edge_weight_median = stats::median(sparse_adj[upper.tri(sparse_adj)][sparse_adj[upper.tri(sparse_adj)] > 0], na.rm = TRUE),
      edge_weight_max = max(sparse_adj[upper.tri(sparse_adj)], na.rm = TRUE),
      graph_edges_nonzero = sum(sparse_adj[upper.tri(sparse_adj)] > 0)
    )
  )
  leiden_backend <- "igraph::cluster_leiden"
} else {
  stop("Unsupported module method in PARAMS$module_method: ", PARAMS$module_method)
}

MEs_train <- moduleEigengenes(train_expr, module_colors)$eigengenes
MEs_train <- orderMEs(MEs_train)

project_to_training_basis <- function(expr_train, expr_new, module_colors, me_train_df) {
  modules <- sort(unique(module_colors))
  valid_scores <- vector("list", length(modules))
  names(valid_scores) <- paste0("ME", modules)
  basis_rows <- vector("list", length(modules))
  for (i in seq_along(modules)) {
    mod <- modules[i]
    me_name <- paste0("ME", mod)
    genes <- names(module_colors)[module_colors == mod]
    x_train <- as.matrix(expr_train[, genes, drop = FALSE])
    x_new <- as.matrix(expr_new[, genes, drop = FALSE])

    mu <- colMeans(x_train, na.rm = TRUE)
    sdv <- apply(x_train, 2, sd, na.rm = TRUE)
    sdv[is.na(sdv) | sdv == 0] <- 1

    z_train <- sweep(sweep(x_train, 2, mu, "-"), 2, sdv, "/")
    z_new <- sweep(sweep(x_new, 2, mu, "-"), 2, sdv, "/")

    if (ncol(z_train) == 1) {
      train_pc1 <- as.numeric(z_train[, 1])
      new_pc1 <- as.numeric(z_new[, 1])
      var_expl <- 1
    } else {
      pca_mod <- prcomp(z_train, center = FALSE, scale. = FALSE)
      load1 <- pca_mod$rotation[, 1]
      train_pc1 <- as.numeric(z_train %*% load1)
      new_pc1 <- as.numeric(z_new %*% load1)
      var_expl <- summary(pca_mod)$importance[2, 1]
    }

    sign_cor <- suppressWarnings(cor(train_pc1, me_train_df[[me_name]], use = "pairwise.complete.obs"))
    sign_flip <- is.finite(sign_cor) && sign_cor < 0
    if (sign_flip) new_pc1 <- -new_pc1

    valid_scores[[me_name]] <- new_pc1
    basis_rows[[i]] <- data.table(
      module = mod,
      eigengene = me_name,
      n_taxa = length(genes),
      pc1_variance_explained = var_expl,
      train_pc1_vs_ME_cor = sign_cor,
      sign_flipped = sign_flip
    )
  }
  list(new_me = as.data.table(valid_scores), basis = rbindlist(basis_rows))
}

MEs_train_dt <- as.data.table(MEs_train)
MEs_train_dt[, sample := rownames(train_expr)]

all_main_ids <- intersect(meta_main$label, rownames(clr))
all_main_expr <- clr[all_main_ids, , drop = FALSE]
proj_all_main <- project_to_training_basis(train_expr, all_main_expr, module_colors, MEs_train)
MEs_all_projected <- proj_all_main$new_me[, names(MEs_train), with = FALSE]
MEs_main_all <- as.data.table(MEs_all_projected)
MEs_main_all[, sample := rownames(all_main_expr)]

valid_expr <- expr_by_core[[PARAMS$validation_core]]
MEs_valid_dt <- MEs_main_all[sample %in% rownames(valid_expr)]

setcolorder(MEs_train_dt, c("sample", setdiff(names(MEs_train_dt), "sample")))
setcolorder(MEs_valid_dt, c("sample", setdiff(names(MEs_valid_dt), "sample")))
setcolorder(MEs_main_all, c("sample", setdiff(names(MEs_main_all), "sample")))

var_r1 <- apply(fit_expr_by_core[["GeoB25202_R1"]], 2, var)
var_r2 <- apply(expr_by_core[[PARAMS$validation_core]], 2, var)
good <- names(which(var_r1 > 0 & var_r2 > 0))

mp <- modulePreservation(
  multiData = list(
    GeoB_R1 = list(data = fit_expr_by_core[["GeoB25202_R1"]][, good, drop = FALSE]),
    GeoB_R2 = list(data = expr_by_core[[PARAMS$validation_core]][, good, drop = FALSE])
  ),
  multiColor = list(GeoB_R1 = module_colors[good]),
  referenceNetworks = 1,
  testNetworks = 2,
  nPermutations = preservation_perms,
  randomSeed = PARAMS$seed,
  verbose = 0
)

pres <- mp$preservation$Z[[1]][[2]]
pres_dt <- data.table(
  module = rownames(pres),
  Zsummary = pres$Zsummary.pres,
  Zdensity = pres$Zdensity.pres,
  Zconnectivity = pres$Zconnectivity.pres
)
pres_dt[, preserved := fcase(Zsummary > 10, "strong", Zsummary > 2, "moderate", default = "weak")]
pres_dt[, module_type := fcase(
  module %in% c("grey", "gold"), "technical",
  default = "biological"
)]
pres_bio_dt <- pres_dt[module_type == "biological"]

age_aligned_concordance <- function(me_all, meta, r1_core, r2_core, validation_core, age_grid_points) {
  me_cols <- grep("^ME", names(me_all), value = TRUE)
  r1_samps <- rownames(fit_expr_by_core[[r1_core]])
  r2_samps <- rownames(expr_by_core[[validation_core]])

  r1_dt <- merge(
    data.table(sample = r1_samps),
    meta[, .(label, y_bp)],
    by.x = "sample",
    by.y = "label",
    all.x = TRUE
  )
  r2_dt <- merge(
    data.table(sample = r2_samps),
    meta[, .(label, y_bp)],
    by.x = "sample",
    by.y = "label",
    all.x = TRUE
  )
  r1_dt[, age_kyr := y_bp / 1000]
  r2_dt[, age_kyr := y_bp / 1000]

  me_r1 <- me_all[sample %in% r1_samps]
  me_r2 <- me_all[sample %in% r2_samps]

  rbindlist(lapply(me_cols, function(me) {
    d1 <- merge(r1_dt[, .(sample, age_kyr)], me_r1[, .(sample, value = get(me))], by = "sample")
    d2 <- merge(r2_dt[, .(sample, age_kyr)], me_r2[, .(sample, value = get(me))], by = "sample")
    d1 <- d1[is.finite(age_kyr) & is.finite(value)][order(age_kyr)]
    d2 <- d2[is.finite(age_kyr) & is.finite(value)][order(age_kyr)]
    if (nrow(d1) < 3 || nrow(d2) < 3) {
      return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    }
    age_min <- max(min(d1$age_kyr), min(d2$age_kyr))
    age_max <- min(max(d1$age_kyr), max(d2$age_kyr))
    if (!is.finite(age_min) || !is.finite(age_max) || age_max <= age_min) {
      return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    }
    xout <- seq(age_min, age_max, length.out = age_grid_points)
    y1 <- approx(d1$age_kyr, d1$value, xout = xout, rule = 2)$y
    y2 <- approx(d2$age_kyr, d2$value, xout = xout, rule = 2)$y
    data.table(
      module = me,
      age_min = age_min,
      age_max = age_max,
      n_grid = age_grid_points,
      pearson_r = suppressWarnings(cor(y1, y2, method = "pearson", use = "complete.obs")),
      spearman_rho = suppressWarnings(cor(y1, y2, method = "spearman", use = "complete.obs")),
      rmse = sqrt(mean((y1 - y2)^2, na.rm = TRUE))
    )
  }), fill = TRUE)
}

concordance_dt <- age_aligned_concordance(
  me_all = MEs_main_all,
  meta = meta_main,
  r1_core = "GeoB25202_R1",
  r2_core = PARAMS$validation_core,
  validation_core = PARAMS$validation_core,
  age_grid_points = PARAMS$wgcna_stability_age_grid_points
)

fwrite(data.table(taxon = names(module_colors), module = module_colors), file.path(DIRS$main, "module_assignments.tsv"), sep = "\t")
fwrite(MEs_main_all, file.path(DIRS$main, "module_eigengenes_main.tsv"), sep = "\t")
fwrite(MEs_train_dt, file.path(DIRS$main, "module_eigengenes_training.tsv"), sep = "\t")
fwrite(MEs_valid_dt, file.path(DIRS$main, "module_eigengenes_validation_projection.tsv"), sep = "\t")
fwrite(proj_all_main$basis, file.path(DIRS$main, "eigengene_projection_basis_main.tsv"), sep = "\t")
fwrite(pres_dt, file.path(DIRS$main, "module_preservation_validation.tsv"), sep = "\t")
fwrite(pres_bio_dt, file.path(DIRS$main, "module_preservation_validation_biological.tsv"), sep = "\t")
fwrite(concordance_dt, file.path(DIRS$main, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
fwrite(sft_dt, file.path(DIRS$main, "soft_power_scan.tsv"), sep = "\t")
saveRDS(net, file.path(DIRS$main, "consensus_wgcna_main.rds"))

module_counts <- data.table(module = unique(module_colors))[order(module)]
module_counts[, n_taxa := as.integer(table(module_colors)[module])]
fwrite(module_counts, file.path(DIRS$main, "module_sizes.tsv"), sep = "\t")

write_run_metadata(
  file.path(DIRS$main, "02_wgcna_main_run_metadata.tsv"),
  "02_wgcna_main.R",
  extra = list(
    training_cores = PARAMS$training_cores,
    validation_core = PARAMS$validation_core,
    wgcna_input_strategy = PARAMS$wgcna_input_strategy,
    wgcna_profile = PARAMS$wgcna_profile,
    wgcna_profile_rationale = PARAMS$wgcna_profile_rationale,
    wgcna_training_samples = nrow(wgcna_training_manifest),
    wgcna_training_counts_by_core = paste(names(fit_core_n), fit_core_n, sep = "=", collapse = ";"),
    hmm_input_balancing_status = PARAMS$hmm_input_balancing_status,
    hmm_input_policy = "all_main_window_samples_projected_from_selected_wgcna_basis",
    module_method = PARAMS$module_method,
    leiden_backend = if (!is.na(leiden_backend)) leiden_backend else "",
    leiden_resolution = PARAMS$leiden_resolution,
    leiden_resolution_strict = if (identical(PARAMS$module_method, "leiden")) PARAMS$leiden_resolution_strict else "",
    leiden_resolution_used = if (identical(PARAMS$module_method, "leiden")) chosen_resolution else "",
    leiden_resolution_auto_adjusted = if (identical(PARAMS$module_method, "leiden")) auto_adjusted else "",
    leiden_graph_source = if (identical(PARAMS$module_method, "leiden")) "WGCNA consensus network" else "",
    leiden_graph_representation = if (identical(PARAMS$module_method, "leiden")) paste0(PARAMS$leiden_graph_symmetrization, " kNN consensus TOM graph") else "",
    leiden_graph_k_neighbors = if (identical(PARAMS$module_method, "leiden")) PARAMS$leiden_graph_k_neighbors else "",
    leiden_graph_symmetrization = if (identical(PARAMS$module_method, "leiden")) PARAMS$leiden_graph_symmetrization else "",
    leiden_consensus_matrix_interpretation = if (identical(PARAMS$module_method, "leiden")) PARAMS$leiden_consensus_matrix_interpretation else "",
    leiden_merge_eigengene_cor = if (identical(PARAMS$module_method, "leiden")) PARAMS$leiden_merge_eigengene_cor else "",
    soft_power = soft_power,
    soft_power_selection = soft_power_selection,
    min_module_size = PARAMS$wgcna_min_module_size,
    deep_split = PARAMS$wgcna_deep_split,
    merge_cut_height = PARAMS$wgcna_merge_cut_height,
    run_mode = run_mode,
    preservation_permutations = preservation_perms,
    stability_age_grid_points = PARAMS$wgcna_stability_age_grid_points,
    results_root = DIRS$results
  )
)
write_session_info(file.path(DIRS$main, "02_wgcna_main_sessionInfo.txt"))

log_msg("02_wgcna_main complete.")
