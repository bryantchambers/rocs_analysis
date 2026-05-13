#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

source(here("config.R"))
set.seed(PARAMS$seed)
options(stringsAsFactors = FALSE)

BASE <- here("networkQC", "input_evaluation")
IN_INPUT <- file.path(BASE, "results", "inputs")
IN_SENS <- file.path(BASE, "results", "sensitivity")
OUT_DEC <- file.path(BASE, "results", "decision")
OUT_FIG <- file.path(BASE, "results", "figures")
REPORT <- file.path(BASE, "INPUT_DECISION_REPORT.md")
dir.create(OUT_DEC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

target_setting <- "exp3"
variants <- c(
  "current_taxon_centered_log",
  "deseq_length_log",
  "sample_clr_raw",
  "log_depth_residualized"
)

log_file <- file.path(OUT_DEC, "progress.log")
if (file.exists(log_file)) invisible(file.remove(log_file))
log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = log_file, append = TRUE)
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

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

state_eta2 <- function(x, state) {
  ok <- is.finite(x) & !is.na(state)
  x <- x[ok]
  state <- as.factor(state[ok])
  if (length(x) < 3 || nlevels(state) < 2) return(NA_real_)
  grand <- mean(x)
  ss_total <- sum((x - grand)^2)
  if (!is.finite(ss_total) || ss_total <= 0) return(NA_real_)
  means <- tapply(x, state, mean)
  sizes <- table(state)
  ss_between <- sum(sizes * (means - grand)^2)
  as.numeric(ss_between / ss_total)
}

clean_variant_name <- function(x) {
  fcase(
    x == "current_taxon_centered_log", "Current exp3",
    x == "deseq_length_log", "DESeq + length log",
    x == "sample_clr_raw", "Sample CLR raw",
    x == "log_depth_residualized", "Depth residualized",
    default = x
  )
}

meta <- fread(file.path(IN_INPUT, "sample_metadata_stage1.tsv"))
depth <- fread(file.path(IN_INPUT, "sample_depth_summary.tsv"))
hmm <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))[, .(sample, state, state_label = label)]

meta_eval <- merge(meta, depth, by.x = "label", by.y = "sample", all.x = TRUE)
meta_eval[, `:=`(
  log_initial = log10(initial + 1),
  log_derep = log10(derep + 1),
  age_kyr = y_bp / 1000
)]
meta_eval <- merge(meta_eval, hmm, by.x = "label", by.y = "sample", all.x = TRUE)

summary_dt <- fread(file.path(IN_SENS, "input_sensitivity_summary.tsv"))[
  setting_id == target_setting & variant %in% variants
]
overlap_dt <- fread(file.path(IN_SENS, "module_overlap_to_current.tsv"))[
  setting_id == target_setting & variant %in% variants
]
pres_dt <- fread(file.path(IN_SENS, "preservation_all.tsv"))[
  setting_id == target_setting & variant %in% variants
]
conc_dt <- fread(file.path(IN_SENS, "eigengene_concordance_all.tsv"))[
  setting_id == target_setting & variant %in% variants
]

current_me <- fread(file.path(IN_SENS, "current_taxon_centered_log", target_setting, "module_eigengenes.tsv"))
current_cols <- setdiff(names(current_me), "sample")
current_module_stats <- NULL

compute_module_stats <- function(me_dt, eval_dt, variant_name) {
  merged <- merge(me_dt, eval_dt, by.x = "sample", by.y = "label", all.x = TRUE)
  me_cols <- setdiff(names(me_dt), "sample")
  rows <- rbindlist(lapply(me_cols, function(me_col) {
    x <- merged[[me_col]]
    data.table(
      variant = variant_name,
      module = sub("^ME", "", me_col),
      cor_d18O = safe_cor(x, merged$mis),
      cor_sst = safe_cor(x, merged$sst),
      cor_log_total = safe_cor(x, merged$log_total_reads),
      cor_library_concentration = safe_cor(x, merged$library_concentration),
      cor_log_initial = safe_cor(x, merged$log_initial),
      cor_log_derep = safe_cor(x, merged$log_derep),
      cor_avg_leng_initial = safe_cor(x, merged$avg_leng_initial),
      cor_avg_len_derep = safe_cor(x, merged$avg_len_derep),
      state_eta2 = state_eta2(x, merged$state_label)
    )
  }), fill = TRUE)
  tech_cols <- grep("^cor_(log_total|library_concentration|log_initial|log_derep|avg_leng_initial|avg_len_derep)$",
                    names(rows), value = TRUE)
  rows[, mean_abs_technical := rowMeans(abs(.SD), na.rm = TRUE), .SDcols = tech_cols]
  rows[, max_abs_technical := do.call(pmax, c(lapply(.SD, abs), na.rm = TRUE)), .SDcols = tech_cols]
  rows
}

log_msg("Computing current exp3 module summaries")
current_module_stats <- compute_module_stats(current_me, meta_eval, "current_taxon_centered_log")
current_module_stats <- current_module_stats[module != "grey"]

module_comp_rows <- list()
tech_rows <- list()
bio_rows <- list()
summary_extra_rows <- list()

for (variant_name in variants) {
  log_msg("Reviewing variant: ", variant_name)
  me_path <- file.path(IN_SENS, variant_name, target_setting, "module_eigengenes.tsv")
  if (!file.exists(me_path)) stop("Missing eigengene file: ", me_path)
  me_dt <- fread(me_path)
  stats_dt <- compute_module_stats(me_dt, meta_eval, variant_name)
  stats_dt <- stats_dt[module != "grey"]

  map_dt <- overlap_dt[variant == variant_name, .(
    current_module = ref_module,
    matched_module = best_test_module,
    overlap_jaccard = best_jaccard
  )]
  map_dt <- map_dt[!is.na(matched_module)]

  comp_rows <- rbindlist(lapply(seq_len(nrow(map_dt)), function(i) {
    row <- map_dt[i]
    current_col <- paste0("ME", row$current_module)
    variant_col <- paste0("ME", row$matched_module)
    if (!current_col %in% names(current_me) || !variant_col %in% names(me_dt)) {
      return(NULL)
    }
    tmp <- merge(
      current_me[, .(sample, current_value = get(current_col))],
      me_dt[, .(sample, variant_value = get(variant_col))],
      by = "sample"
    )
    raw_cor <- safe_cor(tmp$current_value, tmp$variant_value)
    sign_flip <- fifelse(is.finite(raw_cor) & raw_cor < 0, -1, 1)
    aligned_value <- sign_flip * tmp$variant_value

    cur_stat <- current_module_stats[module == row$current_module]
    var_stat <- stats_dt[module == row$matched_module]
    if (nrow(cur_stat) == 0 || nrow(var_stat) == 0) return(NULL)

    data.table(
      variant = variant_name,
      current_module = row$current_module,
      matched_module = row$matched_module,
      overlap_jaccard = row$overlap_jaccard,
      eigengene_corr_raw = raw_cor,
      eigengene_corr_aligned = safe_cor(tmp$current_value, aligned_value),
      sign_flip = sign_flip,
      current_cor_d18O = cur_stat$cor_d18O,
      matched_cor_d18O = sign_flip * var_stat$cor_d18O,
      current_cor_sst = cur_stat$cor_sst,
      matched_cor_sst = sign_flip * var_stat$cor_sst,
      current_state_eta2 = cur_stat$state_eta2,
      matched_state_eta2 = var_stat$state_eta2,
      current_mean_abs_technical = cur_stat$mean_abs_technical,
      matched_mean_abs_technical = var_stat$mean_abs_technical
    )
  }), fill = TRUE)

  if (!nrow(comp_rows)) next

  comp_rows[, `:=`(
    delta_d18O = matched_cor_d18O - current_cor_d18O,
    delta_sst = matched_cor_sst - current_cor_sst,
    delta_state_eta2 = matched_state_eta2 - current_state_eta2,
    delta_mean_abs_technical = matched_mean_abs_technical - current_mean_abs_technical
  )]
  module_comp_rows[[length(module_comp_rows) + 1]] <- comp_rows

  tech_rows[[length(tech_rows) + 1]] <- comp_rows[, .(
    variant,
    current_module,
    matched_module,
    overlap_jaccard,
    eigengene_corr_aligned,
    matched_mean_abs_technical,
    delta_mean_abs_technical
  )]

  bio_rows[[length(bio_rows) + 1]] <- comp_rows[, .(
    variant,
    current_module,
    matched_module,
    overlap_jaccard,
    eigengene_corr_aligned,
    matched_cor_d18O,
    matched_cor_sst,
    matched_state_eta2,
    delta_d18O,
    delta_sst,
    delta_state_eta2,
    sign_flip
  )]

  summary_extra_rows[[length(summary_extra_rows) + 1]] <- comp_rows[, .(
    variant = first(variant),
    mean_aligned_eigengene_corr = mean(eigengene_corr_aligned, na.rm = TRUE),
    min_aligned_eigengene_corr = min(eigengene_corr_aligned, na.rm = TRUE),
    mean_abs_delta_d18O = mean(abs(delta_d18O), na.rm = TRUE),
    mean_abs_delta_sst = mean(abs(delta_sst), na.rm = TRUE),
    mean_abs_delta_state_eta2 = mean(abs(delta_state_eta2), na.rm = TRUE),
    mean_matched_state_eta2 = mean(matched_state_eta2, na.rm = TRUE),
    mean_matched_abs_technical = mean(matched_mean_abs_technical, na.rm = TRUE),
    max_matched_abs_technical = max(matched_mean_abs_technical, na.rm = TRUE)
  )]
}

module_comp_dt <- rbindlist(module_comp_rows, fill = TRUE)
tech_dt <- rbindlist(tech_rows, fill = TRUE)
bio_dt <- rbindlist(bio_rows, fill = TRUE)
summary_extra_dt <- rbindlist(summary_extra_rows, fill = TRUE)

fwrite(module_comp_dt, file.path(OUT_DEC, "module_match_table.tsv"), sep = "\t")
fwrite(tech_dt, file.path(OUT_DEC, "module_technical_correlations.tsv"), sep = "\t")
fwrite(bio_dt, file.path(OUT_DEC, "module_biological_correlations.tsv"), sep = "\t")

pres_summary <- pres_dt[module_type == "biological", .(
  bio_pres_strong_check = sum(preserved == "strong"),
  bio_pres_moderate_check = sum(preserved == "moderate"),
  mean_zsummary = mean(Zsummary, na.rm = TRUE),
  mean_zdensity = mean(Zdensity, na.rm = TRUE),
  mean_zconnectivity = mean(Zconnectivity, na.rm = TRUE)
), by = variant]

conc_summary <- conc_dt[!grepl("grey|gold", module), .(
  mean_pearson_check = mean(pearson_r, na.rm = TRUE),
  mean_spearman_check = mean(spearman_rho, na.rm = TRUE),
  mean_rmse_check = mean(rmse, na.rm = TRUE)
), by = variant]

decision_dt <- merge(summary_dt, summary_extra_dt, by = "variant", all.x = TRUE)
decision_dt <- merge(decision_dt, pres_summary, by = "variant", all.x = TRUE)
decision_dt <- merge(decision_dt, conc_summary, by = "variant", all.x = TRUE)

decision_dt[, score_overlap := norm01(mean_ref_best_jaccard, higher_better = TRUE)]
decision_dt[, score_preservation := rowMeans(cbind(
  norm01(bio_pres_strong, TRUE),
  norm01(mean_concordance_pearson, TRUE),
  norm01(mean_concordance_spearman, TRUE),
  norm01(mean_concordance_rmse, FALSE)
), na.rm = TRUE)]
decision_dt[, score_technical := rowMeans(cbind(
  norm01(abs(cor_PC1_log_total), FALSE),
  norm01(mean_matched_abs_technical, FALSE),
  norm01(max_matched_abs_technical, FALSE)
), na.rm = TRUE)]
decision_dt[, score_biology := rowMeans(cbind(
  norm01(mean_aligned_eigengene_corr, TRUE),
  norm01(mean_abs_delta_d18O, FALSE),
  norm01(mean_abs_delta_sst, FALSE),
  norm01(mean_abs_delta_state_eta2, FALSE),
  norm01(mean_matched_state_eta2, TRUE)
), na.rm = TRUE)]
decision_dt[, score_balance := rowMeans(cbind(
  norm01(grey_pct, FALSE),
  norm01(abs(non_grey_modules - 8), FALSE)
), na.rm = TRUE)]
decision_dt[, decision_score := 0.20 * score_overlap +
              0.20 * score_preservation +
              0.25 * score_technical +
              0.25 * score_biology +
              0.10 * score_balance]

decision_dt[, candidate_class := fcase(
  variant == "log_depth_residualized", "technical-control only",
  variant == "current_taxon_centered_log", "current anchor",
  variant == "deseq_length_log", "candidate alternative",
  variant == "sample_clr_raw", "robustness comparator",
  default = "review"
)]

candidate_dt <- decision_dt[variant != "log_depth_residualized"][order(-decision_score)]
best_candidate <- candidate_dt$variant[[1]]
current_depth_flag <- decision_dt[variant == "current_taxon_centered_log", abs(cor_PC1_log_total)] > 0.8
decision_dt[, recommendation := fcase(
  variant == "log_depth_residualized", "Use only as a technical-control sensitivity.",
  variant == best_candidate & variant == "current_taxon_centered_log" & current_depth_flag,
    "Keep as the operational anchor for now, but the depth signal is still strong; treat `deseq_length_log` as the first fallback if downstream biology looks depth-driven.",
  variant == best_candidate, "Best current candidate among biologically usable inputs.",
  variant == "deseq_length_log", "Closest non-current alternative; review if technical cleanup is needed.",
  variant == "sample_clr_raw", "Important compositional stress test; not the default unless biology is preferred over continuity.",
  variant == "current_taxon_centered_log", "Keep as anchor unless an alternative improves technical behavior without losing biology.",
  default = "Review in context."
)]
decision_dt[, variant_label := clean_variant_name(variant)]
setorder(decision_dt, -decision_score)
decision_dt[, rank := .I]

fwrite(decision_dt, file.path(OUT_DEC, "input_decision_summary.tsv"), sep = "\t")

heat_source <- copy(decision_dt[, .(
  variant_label,
  mean_ref_best_jaccard,
  bio_pres_strong,
  mean_concordance_pearson,
  grey_pct,
  cor_PC1_log_total = abs(cor_PC1_log_total),
  mean_aligned_eigengene_corr,
  mean_abs_delta_d18O,
  mean_abs_delta_sst,
  mean_abs_delta_state_eta2,
  mean_matched_abs_technical
)])
measure_cols <- setdiff(names(heat_source), "variant_label")
heat_source[, (measure_cols) := lapply(.SD, as.numeric), .SDcols = measure_cols]
heat_dt <- melt(
  heat_source,
  id.vars = "variant_label",
  variable.name = "metric",
  value.name = "value"
)
lower_better_metrics <- c(
  "grey_pct", "cor_PC1_log_total", "mean_abs_delta_d18O",
  "mean_abs_delta_sst", "mean_abs_delta_state_eta2", "mean_matched_abs_technical"
)
heat_dt[, fill_value := norm01(value, higher_better = !(first(metric) %in% lower_better_metrics)), by = metric]

p <- ggplot(heat_dt, aes(x = metric, y = variant_label, fill = fill_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient(low = "#b2182b", high = "#2166ac", limits = c(0, 1), guide = "none") +
  labs(
    title = "Input decision metrics for exp3",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(file.path(OUT_FIG, "input_decision_heatmap.png"), p, width = 12, height = 4.5, dpi = 160)

ranked <- decision_dt[order(rank)]
best_row <- ranked[1]
deseq_row <- ranked[variant == "deseq_length_log"]
sample_row <- ranked[variant == "sample_clr_raw"]
resid_row <- ranked[variant == "log_depth_residualized"]

sink(REPORT)
cat("# Input Decision Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: focused decision review for `exp3` across four preprocessing inputs\n\n")

cat("## Plain-language summary\n\n")
cat("The current network uses a taxon-centered log transform. That means each taxon is compared to its own average across samples.\n")
cat("It is useful, but it is not the same thing as a sample-wise CLR, where each sample is centered against all taxa inside that sample.\n\n")
cat("Why this matters: if the transform leaves sample-wide read-depth structure in place, WGCNA may partly organize taxa around technical recovery rather than only around biology.\n\n")
cat("This review asks whether the current `exp3` modules still tell the same biological story when we switch the input frame.\n\n")

cat("## Decision table\n\n")
cat("|rank|input|class|score|mean overlap|grey %|PC1-depth|mean aligned ME corr|mean abs d18O shift|mean abs state shift|\n")
cat("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(ranked))) {
  r <- ranked[i]
  cat(sprintf(
    "|%d|%s|%s|%.3f|%.3f|%.2f|%.3f|%.3f|%.3f|%.3f|\n",
    r$rank, clean_variant_name(r$variant), r$candidate_class, r$decision_score,
    r$mean_ref_best_jaccard, r$grey_pct, abs(r$cor_PC1_log_total),
    r$mean_aligned_eigengene_corr, r$mean_abs_delta_d18O, r$mean_abs_delta_state_eta2
  ))
}
cat("\n")

cat("## Interpretation\n\n")
cat(sprintf("- Best overall usable input in this comparison: `%s`.\n", best_row$variant))
cat(sprintf("- `current exp3`: overlap %.3f, grey %.2f%%, |PC1-depth| %.3f.\n",
            ranked[variant == "current_taxon_centered_log", mean_ref_best_jaccard],
            ranked[variant == "current_taxon_centered_log", grey_pct],
            abs(ranked[variant == "current_taxon_centered_log", cor_PC1_log_total])))
cat(sprintf("- `deseq_length_log`: overlap %.3f, grey %.2f%%, |PC1-depth| %.3f.\n",
            deseq_row$mean_ref_best_jaccard, deseq_row$grey_pct, abs(deseq_row$cor_PC1_log_total)))
cat(sprintf("- `sample_clr_raw`: overlap %.3f, grey %.2f%%, |PC1-depth| %.3f.\n",
            sample_row$mean_ref_best_jaccard, sample_row$grey_pct, abs(sample_row$cor_PC1_log_total)))
cat(sprintf("- `log_depth_residualized`: overlap %.3f, grey %.2f%%, |PC1-depth| %.3f.\n\n",
            resid_row$mean_ref_best_jaccard, resid_row$grey_pct, abs(resid_row$cor_PC1_log_total)))

cat("Read these four inputs this way:\n\n")
cat("- `current exp3`: the continuity anchor. It preserves the existing network by definition and keeps the best overlap.\n")
cat("- `deseq_length_log`: the nearest non-current alternative. If it keeps similar module biology while cleaning technical structure, it is the first real replacement candidate.\n")
cat("- `sample_clr_raw`: the stronger compositional check. If it changes the network a lot, that does not prove it is wrong; it means the choice of compositional frame materially changes the network.\n")
cat("- `log_depth_residualized`: the control. If biology disappears only here, then depth is carrying part of the structure. This is useful diagnostically even if we never use it as the main input.\n\n")

cat("## Recommendation\n\n")
cat(best_row$recommendation, "\n\n")
cat("Do not treat `log_depth_residualized` as the production default unless a later review shows that the current and DESeq-based inputs are clearly dominated by technical structure.\n\n")

cat("## Outputs\n\n")
cat("- Summary: `networkQC/input_evaluation/results/decision/input_decision_summary.tsv`\n")
cat("- Module matching: `networkQC/input_evaluation/results/decision/module_match_table.tsv`\n")
cat("- Technical correlations: `networkQC/input_evaluation/results/decision/module_technical_correlations.tsv`\n")
cat("- Biological correlations: `networkQC/input_evaluation/results/decision/module_biological_correlations.tsv`\n")
cat("- Heatmap: `networkQC/input_evaluation/results/figures/input_decision_heatmap.png`\n")
sink()

log_msg("Input decision review written: ", REPORT)
