#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

OUT_FULL <- here("networkQC", "results", "full_eval")
OUT_FIG <- here("networkQC", "results", "figures")
OUT_TABLE <- here("networkQC", "results", "tables")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(OUT_FULL, "all_settings_summary.tsv")
rank_path <- file.path(OUT_FULL, "all_settings_ranked.tsv")
if (!file.exists(summary_path) || !file.exists(rank_path)) {
  stop("Missing full-eval summary/ranking files.")
}

settings <- fread(summary_path)
ranked <- fread(rank_path)
setorder(ranked, rank)
setting_order <- ranked$setting_id

metric_long <- function(dt, id_cols, metric_cols) {
  melt(dt, id.vars = id_cols, measure.vars = metric_cols,
       variable.name = "metric", value.name = "value")
}

read_setting_modules <- function(setting_id) {
  sdir <- file.path(OUT_FULL, setting_id)
  boot_path <- file.path(sdir, "bootstrap_module_stability_summary.tsv")
  pres_path <- file.path(sdir, "preservation.tsv")
  conc_path <- file.path(sdir, "eigengene_concordance_age_aligned.tsv")
  bal_path <- file.path(sdir, "core_balance_module_jaccard.tsv")

  if (!all(file.exists(c(boot_path, pres_path, conc_path, bal_path)))) {
    return(NULL)
  }

  boot <- fread(boot_path)
  pres <- fread(pres_path)
  conc <- fread(conc_path)
  bal <- fread(bal_path)

  if ("pearson_r.V1" %in% names(conc)) setnames(conc, "pearson_r.V1", "pearson_r")
  if ("spearman_rho.V1" %in% names(conc)) setnames(conc, "spearman_rho.V1", "spearman_rho")
  conc[, module := sub("^ME", "", module)]

  out <- merge(pres, boot, by = "module", all.x = TRUE)
  out <- merge(out, conc, by = "module", all.x = TRUE)
  out <- merge(out, bal, by = "module", all.x = TRUE)
  out[, setting_id := setting_id]
  out[]
}

module_dt <- rbindlist(lapply(settings$setting_id, read_setting_modules), fill = TRUE)
module_dt[, setting_id := factor(setting_id, levels = setting_order)]
module_dt[, module_type := fifelse(is.na(module_type), "unknown", module_type)]
module_dt[, preserved := fifelse(is.na(preserved), "unknown", preserved)]
module_dt[, row_label := sprintf("%s | %s | %s | %s", setting_id, module, module_type, preserved)]

module_metrics <- c(
  "jaccard_median", "jaccard_p05", "jaccard_p95",
  "Zsummary", "Zdensity", "Zconnectivity",
  "pearson_r", "spearman_rho", "rmse", "best_balanced_jaccard"
)
module_heat <- metric_long(module_dt, c("setting_id", "module", "module_type", "preserved", "row_label"), module_metrics)
module_heat[, scaled := {
  r <- range(value, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) rep(0.5, .N) else (value - r[1]) / (r[2] - r[1])
}, by = metric]
module_heat[, label := fifelse(is.na(value), "", sprintf("%.2f", value))]
module_heat[, metric := factor(metric, levels = module_metrics)]
module_heat[, row_label := factor(row_label, levels = rev(unique(module_heat[order(setting_id, module_type, module)]$row_label)))]

fwrite(module_dt, file.path(OUT_TABLE, "full_eval_module_metric_summary.tsv"), sep = "\t")
fwrite(module_heat, file.path(OUT_TABLE, "full_eval_module_metric_heatmap_long.tsv"), sep = "\t")

png(file.path(OUT_FIG, "full_eval_module_metric_heatmap.png"), width = 2600, height = 4200, res = 220)
print(
  ggplot(module_heat, aes(metric, row_label, fill = scaled)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = label), size = 1.55) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey90", name = "Scaled\nwithin metric") +
    labs(
      title = "NetworkQC full-eval module metrics",
      subtitle = "Rows include setting, module, module type, and preservation call; values are printed in original units",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.5),
      plot.title = element_text(face = "bold")
    )
)
dev.off()

setting_metrics <- copy(settings)
boot_agg <- module_dt[module_type == "biological", .(
  median_jaccard_median = median(jaccard_median, na.rm = TRUE),
  median_jaccard_p05 = median(jaccard_p05, na.rm = TRUE),
  median_jaccard_p95 = median(jaccard_p95, na.rm = TRUE),
  median_Zsummary = median(Zsummary, na.rm = TRUE),
  median_Zdensity = median(Zdensity, na.rm = TRUE),
  median_Zconnectivity = median(Zconnectivity, na.rm = TRUE),
  median_module_pearson = median(pearson_r, na.rm = TRUE),
  median_module_spearman = median(spearman_rho, na.rm = TRUE),
  median_module_rmse = median(rmse, na.rm = TRUE),
  median_balanced_jaccard = median(best_balanced_jaccard, na.rm = TRUE),
  weak_bio_modules = sum(preserved == "weak", na.rm = TRUE),
  moderate_bio_modules = sum(preserved == "moderate", na.rm = TRUE),
  strong_bio_modules = sum(preserved == "strong", na.rm = TRUE)
), by = setting_id]

setting_metrics <- merge(setting_metrics, boot_agg, by = "setting_id", all.x = TRUE)
setting_metrics <- merge(setting_metrics, ranked[, .(setting_id, rank, final_score)], by = "setting_id", all.x = TRUE)
setorder(setting_metrics, rank)
fwrite(setting_metrics, file.path(OUT_TABLE, "full_eval_setting_metric_summary.tsv"), sep = "\t")

setting_metric_cols <- c(
  "grey_pct", "non_grey_modules", "bio_pres_strong", "bio_pres_moderate",
  "mean_bootstrap_jaccard", "median_jaccard_p05", "median_jaccard_p95",
  "mean_balanced_jaccard", "mean_concordance_pearson", "mean_concordance_spearman",
  "mean_concordance_rmse", "final_score"
)
setting_heat <- metric_long(setting_metrics, c("setting_id", "rank"), setting_metric_cols)
setting_heat[, scaled := {
  r <- range(value, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) rep(0.5, .N) else (value - r[1]) / (r[2] - r[1])
}, by = metric]
setting_heat[metric %in% c("grey_pct", "mean_concordance_rmse"), scaled := 1 - scaled]
setting_heat[, label := fifelse(is.na(value), "", sprintf("%.2f", value))]
setting_heat[, metric := factor(metric, levels = setting_metric_cols)]
setting_heat[, setting_id := factor(setting_id, levels = rev(setting_order))]
fwrite(setting_heat, file.path(OUT_TABLE, "full_eval_setting_metric_heatmap_long.tsv"), sep = "\t")

png(file.path(OUT_FIG, "full_eval_setting_metric_heatmap.png"), width = 2600, height = 1400, res = 220)
print(
  ggplot(setting_heat, aes(metric, setting_id, fill = scaled)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = 2.1) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f2704", na.value = "grey90", name = "Better\nscaled") +
    labs(
      title = "NetworkQC full-eval setting comparison",
      subtitle = "Lower grey_pct and RMSE are scaled as better; printed values are original units",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold")
    )
)
dev.off()

message("[networkQC] wrote full-eval heatmaps and metric tables")
