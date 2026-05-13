#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
  library(ggplot2)
  library(igraph)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

OUT_BASE <- here("networkQC", "results")
OUT_FULL <- file.path(OUT_BASE, "full_eval")
OUT_TABLE <- file.path(OUT_BASE, "tables")
OUT_FIG <- file.path(OUT_BASE, "figures")
REPORT <- file.path(OUT_BASE, "KME_TOPOLOGY_QC_REPORT.md")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT_BASE, "kme_topology_review.log")
if (file.exists(LOG)) invisible(file.remove(LOG))
log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

norm01 <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  z <- (x - rng[1]) / (rng[2] - rng[1])
  if (!higher_better) z <- 1 - z
  z
}

module_type <- function(x) fifelse(x %in% c("grey", "gold"), "technical", "biological")

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

settings <- data.table(
  setting_id = c("baseline", "opt5", "exp4", "exp3"),
  power = c(20L, 12L, 12L, 12L),
  deepSplit = c(2L, 1L, 3L, 3L),
  mergeCutHeight = c(0.15, 0.20, 0.20, 0.25),
  minModuleSize = c(20L, 20L, 20L, 20L)
)

log_msg("Loading WGCNA input and metadata")
vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
train_expr <- do.call(rbind, lapply(PARAMS$stage1_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
}))
log_msg(sprintf("Training matrix: %d samples x %d taxa", nrow(train_expr), ncol(train_expr)))

read_modules <- function(setting_id) {
  path <- file.path(OUT_FULL, setting_id, "module_assignments.tsv")
  if (!file.exists(path)) stop("Missing module assignment file: ", path)
  mods <- fread(path)
  mods <- mods[taxon %in% colnames(train_expr)]
  setorder(mods, taxon)
  mods
}

compute_kme <- function(setting_id) {
  log_msg("kME pass: ", setting_id)
  mods <- read_modules(setting_id)
  colors <- setNames(mods$module, mods$taxon)
  taxa <- intersect(names(colors), colnames(train_expr))
  expr <- train_expr[, taxa, drop = FALSE]
  colors <- colors[taxa]

  mes <- orderMEs(moduleEigengenes(expr, colors = colors)$eigengenes)
  kme <- cor(expr, mes, use = "pairwise.complete.obs", method = "pearson")
  kme_dt <- as.data.table(kme, keep.rownames = "taxon")
  setnames(kme_dt, old = names(kme_dt), new = sub("^ME", "", names(kme_dt)))

  long <- melt(kme_dt, id.vars = "taxon", variable.name = "module", value.name = "kME")
  long[, module := as.character(module)]
  long[, setting_id := setting_id]

  bio_cols <- setdiff(names(kme_dt), c("taxon", "grey", "gold"))
  kme_dt[, assigned_module := colors[taxon]]
  kme_dt[, assigned_kME := mapply(function(taxon_id, mod) {
    if (!mod %in% names(kme_dt)) return(NA_real_)
    kme_dt[taxon == taxon_id, get(mod)]
  }, taxon, assigned_module)]

  if (length(bio_cols)) {
    kme_dt[, max_bio_kME := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = bio_cols]
    kme_dt[, max_bio_module := bio_cols[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = bio_cols]
    kme_dt[, next_best_bio_kME := mapply(function(taxon_id, assigned) {
      cols <- setdiff(bio_cols, assigned)
      if (!length(cols)) return(NA_real_)
      vals <- unlist(kme_dt[taxon == taxon_id, cols, with = FALSE], use.names = FALSE)
      max(as.numeric(vals), na.rm = TRUE)
    }, taxon, assigned_module)]
  } else {
    kme_dt[, `:=`(max_bio_kME = NA_real_, max_bio_module = NA_character_, next_best_bio_kME = NA_real_)]
  }

  kme_dt[, `:=`(
    setting_id = setting_id,
    module_type = module_type(assigned_module),
    assigned_is_max_kME = assigned_module == max_bio_module,
    kME_margin = assigned_kME - next_best_bio_kME,
    weak_assigned_kME = assigned_kME < 0.2,
    negative_assigned_kME = assigned_kME < 0,
    strong_hub = assigned_kME >= 0.7,
    grey_rescuable = assigned_module == "grey" & max_bio_kME >= 0.5
  )]

  module_summary <- kme_dt[, .(
    n_taxa = .N,
    median_assigned_kME = median(assigned_kME, na.rm = TRUE),
    p05_assigned_kME = safe_quantile(assigned_kME, 0.05),
    frac_assigned_kME_lt_0_2 = mean(weak_assigned_kME, na.rm = TRUE),
    frac_negative_assigned_kME = mean(negative_assigned_kME, na.rm = TRUE),
    frac_assigned_is_max_kME = mean(assigned_is_max_kME, na.rm = TRUE),
    median_kME_margin = median(kME_margin, na.rm = TRUE),
    n_strong_hubs = sum(strong_hub, na.rm = TRUE),
    grey_max_bio_kME_median = if (first(assigned_module) == "grey") median(max_bio_kME, na.rm = TRUE) else NA_real_,
    grey_rescuable_fraction = if (first(assigned_module) == "grey") mean(grey_rescuable, na.rm = TRUE) else NA_real_
  ), by = .(setting_id, module = assigned_module, module_type)]

  setting_summary <- kme_dt[module_type == "biological", .(
    bio_median_assigned_kME = median(assigned_kME, na.rm = TRUE),
    bio_p05_assigned_kME = safe_quantile(assigned_kME, 0.05),
    bio_frac_assigned_kME_lt_0_2 = mean(weak_assigned_kME, na.rm = TRUE),
    bio_frac_negative_assigned_kME = mean(negative_assigned_kME, na.rm = TRUE),
    bio_frac_assigned_is_max_kME = mean(assigned_is_max_kME, na.rm = TRUE),
    bio_median_kME_margin = median(kME_margin, na.rm = TRUE),
    bio_n_strong_hubs = sum(strong_hub, na.rm = TRUE)
  ), by = setting_id]
  grey_summary <- kme_dt[assigned_module == "grey", .(
    grey_max_bio_kME_median = median(max_bio_kME, na.rm = TRUE),
    grey_rescuable_fraction = mean(grey_rescuable, na.rm = TRUE)
  ), by = setting_id]
  setting_summary <- merge(setting_summary, grey_summary, by = "setting_id", all.x = TRUE)

  list(taxon = kme_dt, long = long, module = module_summary, setting = setting_summary)
}

build_tom <- function(power) {
  log_msg("Computing adjacency/TOM for power ", power)
  adj <- adjacency(train_expr, power = power, type = "signed")
  tom <- TOMsimilarity(adj, TOMType = "signed")
  diag(tom) <- 0
  dimnames(tom) <- list(colnames(train_expr), colnames(train_expr))
  log_msg("Finished TOM for power ", power)
  tom
}

threshold_metrics <- function(tom, colors, setting_id, q) {
  taxa <- intersect(names(colors), rownames(tom))
  sub_tom <- tom[taxa, taxa, drop = FALSE]
  vals <- sub_tom[upper.tri(sub_tom)]
  thr <- safe_quantile(vals, q)
  idx <- which(sub_tom >= thr, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  if (!nrow(idx)) {
    return(data.table(
      setting_id = setting_id, q = q, threshold = thr, n_edges = 0L,
      within_edge_fraction = NA_real_, bio_within_edge_fraction = NA_real_,
      modularity = NA_real_, nonisolates = 0L, components = NA_integer_,
      largest_component_fraction = NA_real_, isolate_fraction = 1
    ))
  }

  edges <- data.table(
    from = rownames(sub_tom)[idx[, 1]],
    to = colnames(sub_tom)[idx[, 2]],
    weight = sub_tom[idx]
  )
  edges[, `:=`(
    from_module = colors[from],
    to_module = colors[to]
  )]
  edges[, `:=`(
    from_type = module_type(from_module),
    to_type = module_type(to_module)
  )]
  edges[, `:=`(
    within_module = from_module == to_module,
    bio_edge = from_type == "biological" & to_type == "biological"
  )]

  g <- graph_from_data_frame(edges[, .(from, to, weight)], directed = FALSE,
                             vertices = data.frame(name = taxa))
  V(g)$module <- colors[V(g)$name]
  E(g)$weight <- edges$weight
  deg <- degree(g)
  noniso <- sum(deg > 0)
  comps <- components(g)
  largest <- max(comps$csize)
  membership <- as.integer(factor(V(g)$module))
  mod_val <- modularity(g, membership = membership, weights = E(g)$weight)

  data.table(
    setting_id = setting_id,
    q = q,
    threshold = thr,
    n_edges = nrow(edges),
    within_edge_fraction = mean(edges$within_module, na.rm = TRUE),
    bio_within_edge_fraction = edges[bio_edge == TRUE, mean(within_module, na.rm = TRUE)],
    modularity = mod_val,
    nonisolates = noniso,
    components = comps$no,
    largest_component_fraction = largest / length(taxa),
    isolate_fraction = mean(deg == 0)
  )
}

compute_topology <- function(setting_id, tom) {
  log_msg("Topology pass: ", setting_id)
  mods <- read_modules(setting_id)
  colors <- setNames(mods$module, mods$taxon)
  taxa <- intersect(names(colors), rownames(tom))
  colors <- colors[taxa]
  sub_tom <- tom[taxa, taxa, drop = FALSE]

  bio_taxa <- names(colors)[module_type(colors) == "biological"]
  bio_tom <- sub_tom[bio_taxa, bio_taxa, drop = FALSE]
  pair_idx <- which(upper.tri(bio_tom), arr.ind = TRUE)
  same <- colors[rownames(bio_tom)[pair_idx[, 1]]] == colors[colnames(bio_tom)[pair_idx[, 2]]]
  pair_vals <- bio_tom[pair_idx]
  within_vals <- pair_vals[same]
  between_vals <- pair_vals[!same]
  within_med <- median(within_vals, na.rm = TRUE)
  between_med <- median(between_vals, na.rm = TRUE)

  thresholds <- rbindlist(lapply(c(0.9975, 0.995, 0.99, 0.98), function(q) {
    threshold_metrics(tom, colors, setting_id, q)
  }), fill = TRUE)

  top05 <- thresholds[q == 0.995]
  summary <- data.table(
    setting_id = setting_id,
    within_tom_median = within_med,
    between_tom_median = between_med,
    tom_separation_ratio = within_med / between_med,
    tom_silhouette_like = (within_med - between_med) / max(within_med, between_med, na.rm = TRUE),
    within_edge_fraction_top_0_5pct = top05$within_edge_fraction,
    bio_within_edge_fraction_top_0_5pct = top05$bio_within_edge_fraction,
    modularity_top_0_5pct = top05$modularity,
    largest_component_fraction_top_0_5pct = top05$largest_component_fraction,
    isolate_fraction_top_0_5pct = top05$isolate_fraction
  )
  list(summary = summary, thresholds = thresholds)
}

log_msg("Starting kME review")
kme_results <- lapply(settings$setting_id, compute_kme)
kme_taxon <- rbindlist(lapply(kme_results, `[[`, "taxon"), fill = TRUE)
kme_long <- rbindlist(lapply(kme_results, `[[`, "long"), fill = TRUE)
kme_module <- rbindlist(lapply(kme_results, `[[`, "module"), fill = TRUE)
kme_setting <- rbindlist(lapply(kme_results, `[[`, "setting"), fill = TRUE)

fwrite(kme_taxon, file.path(OUT_TABLE, "kme_taxon_membership.tsv"), sep = "\t")
fwrite(kme_module, file.path(OUT_TABLE, "kme_module_membership_summary.tsv"), sep = "\t")

for (sid in settings$setting_id) {
  dt <- kme_taxon[setting_id == sid & module_type == "biological"]
  p <- ggplot(dt, aes(x = assigned_module, y = assigned_kME, fill = assigned_module)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.75) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "firebrick") +
    coord_cartesian(ylim = c(-0.4, 1)) +
    labs(title = paste("Assigned kME distributions -", sid), x = "Module", y = "Assigned kME") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_FIG, sprintf("kme_module_distributions_%s.png", sid)), p,
         width = 8, height = 4.5, dpi = 160)
}

log_msg("Starting topology review")
tom_by_power <- list()
for (pwr in unique(settings$power)) {
  tom_by_power[[as.character(pwr)]] <- build_tom(pwr)
}
topo_results <- lapply(seq_len(nrow(settings)), function(i) {
  sid <- settings$setting_id[i]
  pwr <- settings$power[i]
  compute_topology(sid, tom_by_power[[as.character(pwr)]])
})
topo_summary <- rbindlist(lapply(topo_results, `[[`, "summary"), fill = TRUE)
topo_thresholds <- rbindlist(lapply(topo_results, `[[`, "thresholds"), fill = TRUE)

fwrite(topo_summary, file.path(OUT_TABLE, "topology_quality_summary.tsv"), sep = "\t")
fwrite(topo_thresholds, file.path(OUT_TABLE, "topology_threshold_sensitivity.tsv"), sep = "\t")

p_thresh <- ggplot(topo_thresholds, aes(x = factor(q), y = within_edge_fraction, color = setting_id, group = setting_id)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.7) +
  scale_x_discrete(labels = c(`0.98` = "top 2%", `0.99` = "top 1%", `0.995` = "top 0.5%", `0.9975` = "top 0.25%")) +
  labs(title = "TOM threshold sensitivity", x = "Edge threshold", y = "Within-module edge fraction") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT_FIG, "topology_threshold_sensitivity.png"), p_thresh, width = 8, height = 4.5, dpi = 160)

log_msg("Building integrated ranking")
full_eval <- fread(file.path(OUT_TABLE, "full_eval_setting_metric_summary.tsv"))
rank_dt <- merge(settings, full_eval[, .(setting_id, final_score, rank, grey_pct, non_grey_modules,
                                         mean_bootstrap_jaccard, mean_balanced_jaccard,
                                         bio_pres_strong, mean_concordance_pearson)],
                 by = "setting_id", all.x = TRUE)
rank_dt <- merge(rank_dt, kme_setting, by = "setting_id", all.x = TRUE)
rank_dt <- merge(rank_dt, topo_summary, by = "setting_id", all.x = TRUE)

rank_dt[, kME_score := rowMeans(cbind(
  norm01(bio_median_assigned_kME, TRUE),
  norm01(bio_p05_assigned_kME, TRUE),
  norm01(bio_frac_assigned_kME_lt_0_2, FALSE),
  norm01(bio_frac_negative_assigned_kME, FALSE),
  norm01(bio_frac_assigned_is_max_kME, TRUE),
  norm01(bio_median_kME_margin, TRUE),
  norm01(grey_rescuable_fraction, FALSE)
), na.rm = TRUE)]

rank_dt[, topology_score := rowMeans(cbind(
  norm01(tom_separation_ratio, TRUE),
  norm01(tom_silhouette_like, TRUE),
  norm01(within_edge_fraction_top_0_5pct, TRUE),
  norm01(bio_within_edge_fraction_top_0_5pct, TRUE),
  norm01(modularity_top_0_5pct, TRUE),
  norm01(isolate_fraction_top_0_5pct, FALSE)
), na.rm = TRUE)]

rank_dt[, integrated_score := 0.55 * norm01(final_score, TRUE) +
          0.25 * kME_score +
          0.20 * topology_score]

rank_dt[, `:=`(
  flag_low_median_kME = bio_median_assigned_kME < 0.50,
  flag_low_p05_kME = bio_p05_assigned_kME < 0.10,
  flag_negative_kME = bio_frac_negative_assigned_kME > 0.05,
  flag_low_assigned_max = bio_frac_assigned_is_max_kME < 0.70,
  flag_low_tom_separation = tom_separation_ratio < 1.25,
  flag_low_within_edges = within_edge_fraction_top_0_5pct < 0.60,
  flag_grey_rescuable = grey_rescuable_fraction > 0.25
)]
flag_cols <- grep("^flag_", names(rank_dt), value = TRUE)
rank_dt[, n_review_flags := rowSums(.SD, na.rm = TRUE), .SDcols = flag_cols]
setorder(rank_dt, -integrated_score)
rank_dt[, integrated_rank := .I]
fwrite(rank_dt, file.path(OUT_TABLE, "final_qc_integrated_ranking.tsv"), sep = "\t")

heat_source <- rank_dt[, .(
  setting_id,
  final_score,
  kME_score,
  topology_score,
  integrated_score,
  grey_pct,
  bio_median_assigned_kME,
  bio_frac_assigned_is_max_kME,
  grey_rescuable_fraction,
  tom_separation_ratio,
  within_edge_fraction_top_0_5pct,
  isolate_fraction_top_0_5pct
)]
measure_cols <- setdiff(names(heat_source), "setting_id")
heat_source[, (measure_cols) := lapply(.SD, as.numeric), .SDcols = measure_cols]
heat_dt <- melt(heat_source, id.vars = "setting_id", variable.name = "metric", value.name = "value")
lower_better <- c("grey_pct", "grey_rescuable_fraction", "isolate_fraction_top_0_5pct")
heat_dt[, fill_value := norm01(value, higher_better = !(first(metric) %in% lower_better)), by = metric]

p_heat <- ggplot(heat_dt, aes(x = metric, y = setting_id, fill = fill_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.8) +
  scale_fill_gradient(low = "#b2182b", high = "#2166ac", limits = c(0, 1), guide = "none") +
  labs(title = "kME and topology QC comparison", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "kme_topology_setting_heatmap.png"), p_heat, width = 12, height = 4.5, dpi = 160)

best <- rank_dt[1]
outcome <- fcase(
  best$setting_id == "exp3", "Keep exp3",
  best$setting_id == "exp4", "Switch to exp4",
  best$setting_id == "opt5", "Hold opt5 as conservative",
  best$setting_id == "baseline", "Baseline recovered unexpectedly; review before accepting",
  default = "Review manually"
)

sink(REPORT)
cat("# kME and Topology QC Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: `baseline`, `opt5`, `exp4`, and `exp3`\n")
cat("- Progress log: `networkQC/results/kme_topology_review.log`\n\n")

cat("## Summary\n\n")
cat("This review adds two WGCNA-native checks to the existing full-eval ranking: kME module membership and TOM topology quality.\n\n")
cat("The input-depth artifact remains a standing caveat. This report chooses the best available module parameters under the current input strategy; it does not claim the depth issue is solved.\n\n")

cat("## Integrated Ranking\n\n")
cat("|rank|setting|full eval|kME|topology|integrated|grey %|bio median kME|TOM separation|flags|\n")
cat("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(rank_dt))) {
  r <- rank_dt[i]
  cat(sprintf("|%d|%s|%.3f|%.3f|%.3f|%.3f|%.2f|%.3f|%.3f|%d|\n",
              r$integrated_rank, r$setting_id, r$final_score, r$kME_score,
              r$topology_score, r$integrated_score, r$grey_pct,
              r$bio_median_assigned_kME, r$tom_separation_ratio, r$n_review_flags))
}
cat("\n")

cat("## Recommendation\n\n")
cat(sprintf("Outcome: `%s`.\n\n", outcome))
cat(sprintf("Top setting: `%s` with integrated score %.3f.\n\n", best$setting_id, best$integrated_score))
cat("Decision rule: prefer `exp3` if it keeps the best combined score without weak kME or topology. Prefer `exp4` only if its kME/topology gain offsets its lower bootstrap/core-balance stability. Treat `opt5` as the conservative fallback if expanded settings look over-fragmented.\n\n")

cat("## Review Flags\n\n")
cat("Flags are warning signs, not automatic rejection rules. They mark settings that deserve inspection before being promoted into the main pipeline.\n\n")
cat("|setting|low median kME|low p05 kME|negative kME|low assigned-is-max|low TOM separation|low within edges|grey rescuable|total flags|\n")
cat("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(rank_dt))) {
  r <- rank_dt[i]
  cat(sprintf("|%s|%s|%s|%s|%s|%s|%s|%s|%d|\n",
              r$setting_id, r$flag_low_median_kME, r$flag_low_p05_kME,
              r$flag_negative_kME, r$flag_low_assigned_max,
              r$flag_low_tom_separation, r$flag_low_within_edges,
              r$flag_grey_rescuable, r$n_review_flags))
}
cat("\n")

cat("## Outputs\n\n")
cat("- kME module summary: `networkQC/results/tables/kme_module_membership_summary.tsv`\n")
cat("- kME taxon table: `networkQC/results/tables/kme_taxon_membership.tsv`\n")
cat("- Topology summary: `networkQC/results/tables/topology_quality_summary.tsv`\n")
cat("- Threshold sensitivity: `networkQC/results/tables/topology_threshold_sensitivity.tsv`\n")
cat("- Integrated ranking: `networkQC/results/tables/final_qc_integrated_ranking.tsv`\n")
cat("- Heatmap: `networkQC/results/figures/kme_topology_setting_heatmap.png`\n")
cat("- Threshold plot: `networkQC/results/figures/topology_threshold_sensitivity.png`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
