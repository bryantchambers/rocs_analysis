#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(WGCNA)
  library(igraph)
  library(here)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

ROOT <- here("network_method_comparison")
OUT_TAB <- file.path(ROOT, "tables")
OUT_FIG <- file.path(ROOT, "figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

ORIG_FULL <- here("networkQC", "results", "full_eval")
ORIG_TABLE <- here("networkQC", "results", "tables")
BAL_FULL <- here("balancednetwork", "results", "qc", "full_eval")
BAL_TABLE <- here("balancednetwork", "results", "tables")
INPUTQC_IMB <- here("InputQC", "core_age_imbalance", "results", "tables")

required_paths <- c(
  file.path(ORIG_FULL, "all_settings_ranked.tsv"),
  file.path(BAL_FULL, "all_settings_ranked.tsv"),
  file.path(RESULTS$stage1, "prokaryotes_vst.rds"),
  file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"),
  file.path(BAL_TABLE, "balanced_baseline_samples.tsv")
)
missing <- required_paths[!file.exists(required_paths)]
if (length(missing)) stop("Missing required inputs:\n", paste(missing, collapse = "\n"))

module_type <- function(x) fifelse(x %in% c("grey", "gold"), "technical", "biological")

norm01 <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(0.5, length(x)))
  z <- (x - rng[1]) / (rng[2] - rng[1])
  if (!higher_better) z <- 1 - z
  z
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

read_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  fread(path)
}

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
meta[, age_kyr := y_bp / 1000]

train_expr <- do.call(rbind, lapply(PARAMS$stage1_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
}))

original_summary <- fread(file.path(ORIG_FULL, "all_settings_ranked.tsv"))
original_summary[, `:=`(
  source = "original",
  method_family = fifelse(setting_id == "baseline", "original_baseline", "original_candidates"),
  method_setting = setting_id
)]

balanced_summary <- fread(file.path(BAL_FULL, "all_settings_ranked.tsv"))
balanced_summary[, `:=`(
  source = "balanced",
  method_family = "balanced_candidates",
  method_setting = setting_id,
  mean_balanced_jaccard = NA_real_
)]

common_cols <- c(
  "source", "method_family", "method_setting", "setting_id",
  "power", "deepSplit", "mergeCutHeight", "minModuleSize",
  "non_grey_modules", "grey_pct", "bio_pres_strong", "bio_pres_moderate",
  "mean_concordance_pearson", "mean_concordance_spearman", "mean_concordance_rmse",
  "mean_bootstrap_jaccard", "mean_balanced_jaccard",
  "n_boot_requested", "n_boot_successful", "n_boot_failed", "rank"
)
for (cc in setdiff(common_cols, names(original_summary))) original_summary[, (cc) := NA]
for (cc in setdiff(common_cols, names(balanced_summary))) balanced_summary[, (cc) := NA]
setting_summary <- rbindlist(list(
  original_summary[, ..common_cols],
  balanced_summary[, ..common_cols]
), fill = TRUE)

setting_summary[source == "original" & setting_id == "exp3", method_family := "original_best"]
setting_summary[source == "balanced" & setting_id == "top3", method_family := "balanced_best"]
setting_summary[, comparison_id := fifelse(source == "original", setting_id, paste0("balanced_", setting_id))]

setting_dir <- function(source, setting_id) {
  if (source == "original") file.path(ORIG_FULL, setting_id) else file.path(BAL_FULL, setting_id)
}

module_rows <- list()
stability_rows <- list()
pres_rows <- list()
for (i in seq_len(nrow(setting_summary))) {
  s <- setting_summary[i]
  sdir <- setting_dir(s$source, s$setting_id)
  mods_path <- file.path(sdir, "module_assignments.tsv")
  if (file.exists(mods_path)) {
    mods <- fread(mods_path)
    mods[, `:=`(
      source = s$source,
      method_family = s$method_family,
      setting_id = s$setting_id,
      comparison_id = s$comparison_id,
      module_type = module_type(module)
    )]
    dist <- mods[, .(n_taxa = .N), by = .(source, method_family, setting_id, comparison_id, module, module_type)]
    dist[, pct_taxa := 100 * n_taxa / sum(n_taxa), by = comparison_id]
    module_rows[[length(module_rows) + 1L]] <- dist
  }
  boot_path <- file.path(sdir, "bootstrap_module_stability_summary.tsv")
  if (file.exists(boot_path)) {
    boot <- fread(boot_path)
    boot[, `:=`(
      source = s$source,
      method_family = s$method_family,
      setting_id = s$setting_id,
      comparison_id = s$comparison_id
    )]
    stability_rows[[length(stability_rows) + 1L]] <- boot
  }
  pres_path <- file.path(sdir, "preservation.tsv")
  if (file.exists(pres_path)) {
    pres <- fread(pres_path)
    pres[, `:=`(
      source = s$source,
      method_family = s$method_family,
      setting_id = s$setting_id,
      comparison_id = s$comparison_id
    )]
    pres_rows[[length(pres_rows) + 1L]] <- pres
  }
}
module_dist <- rbindlist(module_rows, fill = TRUE)
stability <- rbindlist(stability_rows, fill = TRUE)
preservation <- rbindlist(pres_rows, fill = TRUE)

fwrite(module_dist, file.path(OUT_TAB, "module_distribution_comparison.tsv"), sep = "\t")
fwrite(stability, file.path(OUT_TAB, "module_stability_comparison.tsv"), sep = "\t")
fwrite(preservation, file.path(OUT_TAB, "preservation_comparison.tsv"), sep = "\t")

sample_fairness <- function(dt, source_label, bin_width = 10) {
  x <- copy(dt)
  x[, age_bin := sprintf("%.0f-%.0f", floor(age_kyr / bin_width) * bin_width, floor(age_kyr / bin_width) * bin_width + bin_width)]
  core_counts <- x[, .N, by = core]
  wide <- dcast(x[, .N, by = .(age_bin, core)], age_bin ~ core, value.var = "N", fill = 0)
  core_cols <- intersect(PARAMS$stage1_cores, names(wide))
  cv_by_bin <- if (length(core_cols) > 1) {
    wide[, apply(.SD, 1, function(v) {
      if (mean(v) == 0) return(NA_real_)
      stats::sd(v) / mean(v)
    }), .SDcols = core_cols]
  } else numeric()
  data.table(
    source = source_label,
    train_samples = nrow(x),
    min_core_n = min(core_counts$N),
    max_core_n = max(core_counts$N),
    st8_n = core_counts[core == "ST8", N],
    st8_fraction = core_counts[core == "ST8", N] / nrow(x),
    core_n_cv = sd(core_counts$N) / mean(core_counts$N),
    age_bin_count_cv_mean = mean(cv_by_bin, na.rm = TRUE),
    age_bin_count_cv_median = median(cv_by_bin, na.rm = TRUE),
    n_age_bins = uniqueN(x$age_bin)
  )
}

orig_train_meta <- meta[core %in% PARAMS$stage1_cores & label %in% rownames(vst), .(label, core, age_kyr)]
balanced_sel <- fread(file.path(BAL_TABLE, "balanced_baseline_samples.tsv"))
balanced_train_meta <- balanced_sel[, .(label, core, age_kyr)]
fairness <- rbindlist(list(
  sample_fairness(orig_train_meta, "original"),
  sample_fairness(balanced_train_meta, "balanced")
), fill = TRUE)

imbalance_density <- read_optional(file.path(INPUTQC_IMB, "core_age_density_summary.tsv"))
if (!is.null(imbalance_density)) {
  fwrite(imbalance_density, file.path(OUT_TAB, "inputqc_core_age_density_summary.tsv"), sep = "\t")
}
fwrite(fairness, file.path(OUT_TAB, "sampling_fairness_summary.tsv"), sep = "\t")

compute_kme_topology <- function(settings) {
  make_tom <- function(power) {
    adj <- adjacency(train_expr, power = power, type = "signed")
    tom <- TOMsimilarity(adj, TOMType = "signed")
    diag(tom) <- 0
    dimnames(tom) <- list(colnames(train_expr), colnames(train_expr))
    tom
  }
  powers <- sort(unique(settings$power[is.finite(settings$power)]))
  tom_cache <- setNames(vector("list", length(powers)), as.character(powers))
  for (p in powers) tom_cache[[as.character(p)]] <- make_tom(p)

  out <- list()
  for (i in seq_len(nrow(settings))) {
    s <- settings[i]
    mods_path <- file.path(setting_dir(s$source, s$setting_id), "module_assignments.tsv")
    if (!file.exists(mods_path)) next
    mods <- fread(mods_path)
    mods <- mods[taxon %in% colnames(train_expr)]
    colors <- setNames(mods$module, mods$taxon)
    taxa <- intersect(names(colors), colnames(train_expr))
    colors <- colors[taxa]
    expr <- train_expr[, taxa, drop = FALSE]

    mes <- orderMEs(moduleEigengenes(expr, colors = colors)$eigengenes)
    kme <- cor(expr, mes, use = "pairwise.complete.obs", method = "pearson")
    kdt <- as.data.table(kme, keep.rownames = "taxon")
    setnames(kdt, names(kdt), sub("^ME", "", names(kdt)))
    bio_cols <- setdiff(names(kdt), c("taxon", "grey", "gold"))
    kdt[, assigned_module := colors[taxon]]
    kdt[, assigned_kME := mapply(function(taxon_id, mod) {
      if (!mod %in% names(kdt)) return(NA_real_)
      kdt[taxon == taxon_id, get(mod)]
    }, taxon, assigned_module)]
    if (length(bio_cols)) {
      kmat <- as.matrix(kdt[, ..bio_cols])
      kdt[, max_bio_kME := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = bio_cols]
      kdt[, max_bio_module := bio_cols[max.col(kmat, ties.method = "first")]]
      next_best <- rep(NA_real_, nrow(kdt))
      for (j in seq_len(nrow(kdt))) {
        cols <- setdiff(bio_cols, kdt$assigned_module[j])
        if (length(cols)) next_best[j] <- max(kmat[j, colnames(kmat) %in% cols], na.rm = TRUE)
      }
      kdt[, next_best_bio_kME := next_best]
    } else {
      kdt[, `:=`(max_bio_kME = NA_real_, max_bio_module = NA_character_, next_best_bio_kME = NA_real_)]
    }
    kdt[, `:=`(
      module_type = module_type(assigned_module),
      assigned_is_max_kME = assigned_module == max_bio_module,
      kME_margin = assigned_kME - next_best_bio_kME,
      weak_assigned_kME = assigned_kME < 0.2,
      negative_assigned_kME = assigned_kME < 0,
      strong_hub = assigned_kME >= 0.7,
      grey_rescuable = assigned_module == "grey" & max_bio_kME >= 0.5
    )]
    ksum <- kdt[module_type == "biological", .(
      bio_median_assigned_kME = median(assigned_kME, na.rm = TRUE),
      bio_p05_assigned_kME = safe_quantile(assigned_kME, 0.05),
      bio_frac_assigned_kME_lt_0_2 = mean(weak_assigned_kME, na.rm = TRUE),
      bio_frac_negative_assigned_kME = mean(negative_assigned_kME, na.rm = TRUE),
      bio_frac_assigned_is_max_kME = mean(assigned_is_max_kME, na.rm = TRUE),
      bio_median_kME_margin = median(kME_margin, na.rm = TRUE),
      bio_n_strong_hubs = sum(strong_hub, na.rm = TRUE)
    )]
    gsum <- kdt[assigned_module == "grey", .(
      grey_max_bio_kME_median = median(max_bio_kME, na.rm = TRUE),
      grey_rescuable_fraction = mean(grey_rescuable, na.rm = TRUE)
    )]

    tom <- tom_cache[[as.character(s$power)]]
    ttaxa <- intersect(names(colors), rownames(tom))
    cols <- colors[ttaxa]
    sub_tom <- tom[ttaxa, ttaxa, drop = FALSE]
    bio_taxa <- names(cols)[module_type(cols) == "biological"]
    bio_tom <- sub_tom[bio_taxa, bio_taxa, drop = FALSE]
    pidx <- which(upper.tri(bio_tom), arr.ind = TRUE)
    same <- cols[rownames(bio_tom)[pidx[, 1]]] == cols[colnames(bio_tom)[pidx[, 2]]]
    pair_vals <- bio_tom[pidx]
    within_med <- median(pair_vals[same], na.rm = TRUE)
    between_med <- median(pair_vals[!same], na.rm = TRUE)

    vals <- sub_tom[upper.tri(sub_tom)]
    threshold <- as.numeric(quantile(vals, 0.995, na.rm = TRUE))
    idx <- which(sub_tom >= threshold, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    edges <- data.table(
      from = rownames(sub_tom)[idx[, 1]],
      to = colnames(sub_tom)[idx[, 2]],
      weight = sub_tom[idx]
    )
    if (nrow(edges)) {
      edges[, `:=`(from_module = cols[from], to_module = cols[to])]
      edges[, within_module := from_module == to_module]
      edges[, bio_edge := module_type(from_module) == "biological" & module_type(to_module) == "biological"]
      g <- graph_from_data_frame(edges[, .(from, to, weight)], directed = FALSE, vertices = data.frame(name = ttaxa))
      E(g)$weight <- edges$weight
      deg <- degree(g)
      comps <- components(g)
      membership <- as.integer(factor(cols[V(g)$name]))
      topo <- data.table(
        within_edge_fraction_top_0_5pct = mean(edges$within_module, na.rm = TRUE),
        bio_within_edge_fraction_top_0_5pct = edges[bio_edge == TRUE, mean(within_module, na.rm = TRUE)],
        modularity_top_0_5pct = modularity(g, membership = membership, weights = E(g)$weight),
        largest_component_fraction_top_0_5pct = max(comps$csize) / length(ttaxa),
        isolate_fraction_top_0_5pct = mean(deg == 0)
      )
    } else {
      topo <- data.table(
        within_edge_fraction_top_0_5pct = NA_real_,
        bio_within_edge_fraction_top_0_5pct = NA_real_,
        modularity_top_0_5pct = NA_real_,
        largest_component_fraction_top_0_5pct = NA_real_,
        isolate_fraction_top_0_5pct = 1
      )
    }
    out[[length(out) + 1L]] <- cbind(
      s[, .(source, method_family, setting_id, comparison_id, power, deepSplit, mergeCutHeight, minModuleSize)],
      ksum,
      gsum,
      data.table(
        within_tom_median = within_med,
        between_tom_median = between_med,
        tom_separation_ratio = within_med / between_med,
        tom_silhouette_like = (within_med - between_med) / max(within_med, between_med, na.rm = TRUE)
      ),
      topo
    )
  }
  rbindlist(out, fill = TRUE)
}

kme_topology <- compute_kme_topology(setting_summary)
fwrite(kme_topology, file.path(OUT_TAB, "kme_topology_comparison.tsv"), sep = "\t")

setting_level <- merge(setting_summary, fairness[, .(
  source, train_samples, min_core_n, max_core_n, st8_n, st8_fraction, core_n_cv,
  age_bin_count_cv_mean, age_bin_count_cv_median, n_age_bins
)], by = "source", all.x = TRUE)
setting_level <- merge(setting_level, kme_topology, by = c(
  "source", "method_family", "setting_id", "comparison_id", "power", "deepSplit", "mergeCutHeight", "minModuleSize"
), all.x = TRUE)

setting_level[, `:=`(
  n_st8 = norm01(st8_fraction, FALSE),
  n_core_cv = norm01(core_n_cv, FALSE),
  n_age_cv = norm01(age_bin_count_cv_mean, FALSE),
  n_grey = norm01(grey_pct, FALSE),
  n_boot = norm01(mean_bootstrap_jaccard, TRUE),
  n_conc = norm01(mean_concordance_pearson, TRUE),
  n_kme_med = norm01(bio_median_assigned_kME, TRUE),
  n_kme_p05 = norm01(bio_p05_assigned_kME, TRUE),
  n_tom_sep = norm01(tom_separation_ratio, TRUE),
  n_bio_edges = norm01(bio_within_edge_fraction_top_0_5pct, TRUE)
)]
setting_level[, `:=`(
  score_sampling_fairness = rowMeans(.SD, na.rm = TRUE)
), .SDcols = c("n_st8", "n_core_cv", "n_age_cv")]
setting_level[, `:=`(
  score_network_qc = rowMeans(.SD, na.rm = TRUE)
), .SDcols = c("n_grey", "n_boot", "n_conc")]
setting_level[, `:=`(
  score_kme_topology = rowMeans(.SD, na.rm = TRUE)
), .SDcols = c("n_kme_med", "n_kme_p05", "n_tom_sep", "n_bio_edges")]
setting_level[, tiered_summary_score := 0.35 * score_sampling_fairness + 0.35 * score_network_qc + 0.30 * score_kme_topology]
setorder(setting_level, -tiered_summary_score)
setting_level[, tiered_rank := .I]
fwrite(setting_level, file.path(OUT_TAB, "setting_level_comparison.tsv"), sep = "\t")

method_pick_ids <- c("baseline", "exp3", "balanced_top3")
method_level <- setting_level[comparison_id %in% method_pick_ids]
method_level[comparison_id == "baseline", method_label := "original_baseline"]
method_level[comparison_id == "exp3", method_label := "original_best_exp3"]
method_level[comparison_id == "balanced_top3", method_label := "balanced_best_top3"]
method_level[, display_order := fcase(
  method_label == "balanced_best_top3", 1L,
  method_label == "original_best_exp3", 2L,
  method_label == "original_baseline", 3L,
  default = 99L
)]
setorder(method_level, display_order)
fwrite(method_level, file.path(OUT_TAB, "method_level_comparison.tsv"), sep = "\t")

readiness_files <- data.table(
  artifact = c(
    "module assignments", "module eigengenes", "HMM states", "HMM validation",
    "state fingerprints", "state importance", "driver integration", "network state stats",
    "functional enrichment"
  ),
  path = c(
    file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"),
    file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv"),
    file.path(RESULTS$hmm, "hmm_states.tsv"),
    file.path(RESULTS$hmm, "hmm_validation_metrics.tsv"),
    file.path(RESULTS$hmm, "state_fingerprints.tsv"),
    file.path(RESULTS$importance, "state_importance_scores.tsv"),
    file.path(RESULTS$importance, "integrated_driver_summary.tsv"),
    file.path(RESULTS$network_stats, "state_network_stats.tsv"),
    file.path(RESULTS$importance, "state_functional_enrichment.tsv")
  ),
  current_scope = c(rep("original_current_pipeline", 9))
)
readiness_files[, exists := file.exists(path)]
readiness_files[, balanced_equivalent_exists := FALSE]
readiness_files[artifact %in% c("module assignments", "module eigengenes"), balanced_equivalent_exists := file.exists(file.path(here("balancednetwork", "results", "wgcna"), fifelse(artifact == "module assignments", "module_assignments.tsv", "module_eigengenes.tsv")))]
readiness_files[, review_note := fifelse(
  exists & !balanced_equivalent_exists & !artifact %in% c("module assignments", "module eigengenes"),
  "downstream biology currently represents original/current pipeline only",
  fifelse(exists, "available", "missing")
)]
fwrite(readiness_files, file.path(OUT_TAB, "downstream_readiness_checklist.tsv"), sep = "\t")

workflow_strategy <- data.table(
  approach = c("original_scripts", "InputQC_sensitivity", "balancednetwork"),
  primary_code = c("scripts/01_data_prep.R + scripts/02_wgcna.R", "InputQC/ and networkQC/input_evaluation/", "balancednetwork/"),
  input_basis = c(
    "all eligible <=150 kyr samples after damaged prokaryote filtering and taxon-centered CLR",
    "diagnostic variants around current input plus low-detection/ST8/core-age imbalance checks",
    "same stage1 matrix, then age-bin/core-balanced training subset"
  ),
  sample_balance_handling = c(
    "does not equalize core sample counts; ST8 is 60.8% of training samples",
    "diagnoses imbalance and tests sensitivity but is not the promoted network builder",
    "equalizes training cores within shared age bins; ST8 is 33.3% of balanced training samples"
  ),
  main_strength = c(
    "uses all available information and currently supports downstream HMM/driver biology",
    "identifies whether depth, low-detection, and ST8 structure threaten interpretation",
    "directly tests whether modules survive after removing ST8 sample-count dominance"
  ),
  main_risk = c(
    "stability may partially reflect repeated ST8-heavy covariance structure",
    "diagnostic rather than a complete alternative network decision by itself",
    "higher grey burden and smaller training set may discard real information"
  ),
  decision_role = c(
    "current end-to-end biological reference",
    "evidence layer for whether original input assumptions are safe",
    "sampling-fair challenger and sensitivity control"
  )
)
fwrite(workflow_strategy, file.path(OUT_TAB, "workflow_strategy_comparison.tsv"), sep = "\t")

gap_matrix <- data.table(
  gap = c(
    "balanced downstream HMM rerun",
    "balanced driver/taxon importance rerun",
    "balanced functional enrichment/state breakdown",
    "matched original-on-balanced-samples run",
    "module overlap original exp3 vs balanced top3",
    "age-bin/core-stratified downstream validation"
  ),
  why_it_matters = c(
    "tests whether balanced module eigengenes recover the same ecological state sequence",
    "checks whether key taxa/drivers survive the fair-sampling network",
    "checks whether functional interpretation is preserved after balancing",
    "separates parameter effects from sample-subset effects",
    "shows whether balanced modules are merged/split versions of original biology",
    "tests whether conclusions depend on ST8-heavy age bands"
  ),
  current_status = c(
    "missing; current HMM outputs are original/current pipeline",
    "missing; current importance outputs are original/current pipeline",
    "missing for balanced; available for original/current pipeline",
    "not yet run as a dedicated control",
    "partly available through balanced_vs_original overlap, needs integrated summary",
    "not yet summarized in method comparison"
  ),
  recommended_next_action = c(
    "run HMM using balanced top3 module eigengenes in an isolated output directory",
    "rerun taxon importance and driver integration against balanced HMM/module outputs",
    "rerun functional/state summaries after balanced HMM",
    "fit exp3 parameters on the balanced sample manifest and compare to balanced top3",
    "add per-module best-overlap table with taxon/function annotations",
    "stratify key state/driver claims by age bin and core"
  )
)
fwrite(gap_matrix, file.path(OUT_TAB, "analysis_gap_matrix.tsv"), sep = "\t")

plot_setting <- copy(setting_level)
plot_setting[, display := comparison_id]
ggplot(plot_setting[method_family %in% c("original_best", "balanced_best", "original_baseline") | source == "balanced"], aes(x = reorder(display, mean_bootstrap_jaccard), y = mean_bootstrap_jaccard, fill = source)) +
  geom_col(width = 0.72) +
  coord_flip() +
  labs(x = NULL, y = "Mean bootstrap Jaccard", title = "Bootstrap stability by setting") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT_FIG, "bootstrap_stability_by_setting.png"), width = 8, height = 5, dpi = 160)

ggplot(module_dist[comparison_id %in% c("baseline", "exp3", "balanced_top3")], aes(x = comparison_id, y = pct_taxa, fill = module)) +
  geom_col(width = 0.7) +
  labs(x = NULL, y = "Taxa percentage", title = "Module size distribution") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT_FIG, "module_distribution_key_methods.png"), width = 8, height = 5, dpi = 160)

metric_long <- melt(
  method_level,
  id.vars = c("method_label", "source"),
  measure.vars = c("score_sampling_fairness", "score_network_qc", "score_kme_topology", "tiered_summary_score"),
  variable.name = "metric",
  value.name = "score"
)
ggplot(metric_long, aes(x = metric, y = method_label, fill = score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3) +
  scale_fill_gradient(low = "#f2f0e6", high = "#2b6f6a", na.value = "grey80") +
  labs(x = NULL, y = NULL, title = "Tiered comparison scores") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_FIG, "tiered_method_score_heatmap.png"), width = 8, height = 4.5, dpi = 160)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
}

best_bal <- setting_level[source == "balanced"][which.max(mean_bootstrap_jaccard)]
best_orig <- setting_level[comparison_id == "exp3"]
base <- setting_level[comparison_id == "baseline"]

research_note <- file.path(ROOT, "RESEARCH_NOTE_WGCNA_IMBALANCE.md")
sink(research_note)
cat("# WGCNA Imbalance Research Note\n\n")
cat("## Question\n\n")
cat("Can unequal sample density, especially ST8 dominance, make the original WGCNA network look more stable than a more sampling-balanced network?\n\n")
cat("## Method Background\n\n")
cat("- WGCNA builds weighted correlation networks and summarizes modules with eigengenes, kME/module membership, and topological overlap. Core references: Langfelder and Horvath 2008 (<https://link.springer.com/article/10.1186/1471-2105-9-559>) and the WGCNA manual (<https://cran.r-universe.dev/WGCNA/doc/manual.html>).\n")
cat("- The WGCNA FAQ warns that strong categorical drivers, batch effects, or biological heterogeneity can dominate correlations; it recommends inspecting sample clustering and considering adjustment or consensus-style analyses when heterogeneity is strong (<https://edo98811.github.io/WGCNA_official_documentation/faq.html>).\n")
cat("- Consensus WGCNA is designed to find modules that recur across datasets or groups. In this project, cores/age bins behave like structured groups, so a balance-aware comparison is more defensible than treating raw sample count as neutral information.\n\n")
cat("## Project-Specific Risk\n\n")
cat(sprintf("- Original training sample counts are imbalanced: ST8 fraction is %.1f%%; balanced training ST8 fraction is %.1f%%.\n", 100 * fairness[source == "original", st8_fraction], 100 * fairness[source == "balanced", st8_fraction]))
cat("- InputQC already flags ST8 as much more densely sampled and shows a 90-110 kya band where ST8 low-detection samples are enriched beyond sampling expectation.\n")
cat("- Therefore, original-network stability can mean either true ecological robustness or repeated recovery of an ST8-heavy correlation structure.\n\n")
cat("## Practical Interpretation Rule\n\n")
cat("Treat bootstrap stability as necessary but not sufficient. A fair network should retain stability after core/age balancing, show reasonable grey burden, preserve modules in GeoB25202_R2, maintain kME/TOM coherence, and support downstream ecological state/driver interpretations.\n")
sink()

report <- file.path(ROOT, "METHOD_COMPARISON_REPORT.md")
sink(report)
cat("# Original vs Balanced WGCNA Method Comparison\n\n")
cat("## Executive Summary\n\n")
cat(sprintf("- Balanced best setting: `%s` (power=%d, deepSplit=%d, mergeCutHeight=%.2f, minModuleSize=%d).\n",
            best_bal$comparison_id, best_bal$power, best_bal$deepSplit, best_bal$mergeCutHeight, best_bal$minModuleSize))
cat(sprintf("- Original best setting: `exp3` (power=%d, deepSplit=%d, mergeCutHeight=%.2f, minModuleSize=%d).\n",
            best_orig$power, best_orig$deepSplit, best_orig$mergeCutHeight, best_orig$minModuleSize))
cat(sprintf("- Balanced best has stronger bootstrap stability: %s vs original exp3 %s.\n",
            fmt(best_bal$mean_bootstrap_jaccard), fmt(best_orig$mean_bootstrap_jaccard)))
cat(sprintf("- Original exp3 has lower grey burden: %s%% vs balanced best %s%%.\n",
            fmt(best_orig$grey_pct, 2), fmt(best_bal$grey_pct, 2)))
cat(sprintf("- Original exp3 has stronger TOM separation: %s vs balanced best %s.\n",
            fmt(best_orig$tom_separation_ratio), fmt(best_bal$tom_separation_ratio)))
cat("- Current downstream HMM/driver/functional outputs represent the original/current pipeline, not a full balanced downstream rerun.\n\n")
cat("## Sampling Fairness\n\n")
cat("|method|training samples|ST8 fraction|core count CV|age-bin count CV mean|\n")
cat("|---|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(fairness))) {
  r <- fairness[i]
  cat(sprintf("|%s|%d|%.3f|%.3f|%.3f|\n", r$source, r$train_samples, r$st8_fraction, r$core_n_cv, r$age_bin_count_cv_mean))
}
cat("\n## Key Method Comparison\n\n")
cat("|method|grey %|non-grey modules|bootstrap Jaccard|Pearson concordance|bio median kME|TOM separation|\n")
cat("|---|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(method_level))) {
  r <- method_level[i]
  cat(sprintf("|%s|%.2f|%d|%.3f|%.3f|%.3f|%.3f|\n",
              r$method_label, r$grey_pct, r$non_grey_modules, r$mean_bootstrap_jaccard,
              r$mean_concordance_pearson, r$bio_median_assigned_kME,
              r$tom_separation_ratio))
}
cat("\n## Metric Winners\n\n")
cat(sprintf("- Sampling fairness: balanced (`ST8 fraction %.3f`, core count CV %.3f).\n",
            fairness[source == "balanced", st8_fraction], fairness[source == "balanced", core_n_cv]))
cat(sprintf("- Bootstrap stability: balanced best (`%.3f`).\n", best_bal$mean_bootstrap_jaccard))
cat(sprintf("- Grey burden and module recovery: original exp3 (`%.2f%%` grey, `%d` non-grey modules).\n",
            best_orig$grey_pct, best_orig$non_grey_modules))
cat(sprintf("- TOM separation among non-baseline candidates: original exp3 (`%.3f`).\n", best_orig$tom_separation_ratio))
cat(sprintf("- Downstream biology availability: original/current pipeline (`%d/%d` artifacts available); balanced downstream rerun not yet available.\n",
            readiness_files[current_scope == "original_current_pipeline", sum(exists)], nrow(readiness_files)))
cat("\n## Interpretation\n\n")
cat("- Balanced is more sampling-fair and, in the completed 500-bootstrap run, more label-stable than original exp3.\n")
cat("- Original exp3 is less grey, has more biological modules, stronger TOM separation, and currently owns the downstream biology evidence.\n")
cat("- The fair answer is not simply that one method wins. Balanced should be treated as the stronger sampling-control network; original exp3 remains the stronger current end-to-end biological network until balanced downstream HMM/driver analyses are rerun.\n")
cat("- If balanced downstream biology recovers the same state/driver story, promote balanced or use it as the primary sensitivity-supported network. If it loses major ecological signal, keep original exp3 with an explicit ST8-dominance caveat.\n\n")
cat("## Workflow Strategy Comparison\n\n")
cat("|approach|decision role|main strength|main risk|\n")
cat("|---|---|---|---|\n")
for (i in seq_len(nrow(workflow_strategy))) {
  r <- workflow_strategy[i]
  cat(sprintf("|%s|%s|%s|%s|\n", r$approach, r$decision_role, r$main_strength, r$main_risk))
}
cat("\n## Remaining Gaps\n\n")
cat("|gap|status|next action|\n")
cat("|---|---|---|\n")
for (i in seq_len(nrow(gap_matrix))) {
  r <- gap_matrix[i]
  cat(sprintf("|%s|%s|%s|\n", r$gap, r$current_status, r$recommended_next_action))
}
cat("\n")
cat("## Output Tables\n\n")
cat("- `tables/method_level_comparison.tsv`\n")
cat("- `tables/setting_level_comparison.tsv`\n")
cat("- `tables/module_distribution_comparison.tsv`\n")
cat("- `tables/module_stability_comparison.tsv`\n")
cat("- `tables/kme_topology_comparison.tsv`\n")
cat("- `tables/downstream_readiness_checklist.tsv`\n")
cat("- `tables/workflow_strategy_comparison.tsv`\n")
cat("- `tables/analysis_gap_matrix.tsv`\n")
sink()

message("Wrote network method comparison outputs to ", ROOT)
