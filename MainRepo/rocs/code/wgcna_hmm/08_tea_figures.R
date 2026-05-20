#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)

INDEX_COLS <- c("OAP", "MGI", "MFI", "DCI", "MII", "SRPI")
OUT_DIR <- file.path(DIRS$reports, "tea_figures")
OUT_TBL <- file.path(OUT_DIR, "tables")

ensure_dirs(c(OUT_DIR, OUT_TBL))

required_inputs <- c(
  file.path(DIRS$main_tea, "tea_indices_per_sample_module_main.tsv"),
  file.path(DIRS$main_tea, "tea_indices_per_sample_main.tsv"),
  file.path(DIRS$main_tea, "tea_module_centered_values.tsv"),
  file.path(DIRS$main_tea, "tea_state_summary_main.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_nondamaged_indices_per_sample_main.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_temporal_descriptive_summary.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_state_module_variability_specialization_summary.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_dominant_metabolism_calls.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_nondamaged_vs_module_comparison.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_nondamaged_vs_module_by_index.tsv"),
  file.path(DIRS$main, "tea_variability", "tea_brown_vs_nondamaged_comparison.tsv"),
  file.path(DIRS$main, "sample_metadata_main.tsv"),
  file.path(DIRS$main, "hmm_states_main.tsv"),
  file.path(DIRS$main, "state_fingerprints_main.tsv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required inputs for 08_tea_figures: ", paste(missing_inputs, collapse = "; "))
}

log_msg("Loading TEA inputs for exploratory figure battery...")

tea_sample_module <- fread(file.path(DIRS$main_tea, "tea_indices_per_sample_module_main.tsv"))
tea_sample <- fread(file.path(DIRS$main_tea, "tea_indices_per_sample_main.tsv"))
tea_centered <- fread(file.path(DIRS$main_tea, "tea_module_centered_values.tsv"))
tea_state_summary <- fread(file.path(DIRS$main_tea, "tea_state_summary_main.tsv"))
nondmg_sample <- fread(file.path(DIRS$main, "tea_variability", "tea_nondamaged_indices_per_sample_main.tsv"))
temporal_summary <- fread(file.path(DIRS$main, "tea_variability", "tea_temporal_descriptive_summary.tsv"))
state_module_specialization <- fread(file.path(DIRS$main, "tea_variability", "tea_state_module_variability_specialization_summary.tsv"))
dominance_calls <- fread(file.path(DIRS$main, "tea_variability", "tea_dominant_metabolism_calls.tsv"))
nondmg_vs_module <- fread(file.path(DIRS$main, "tea_variability", "tea_nondamaged_vs_module_comparison.tsv"))
nondmg_vs_module_by_index <- fread(file.path(DIRS$main, "tea_variability", "tea_nondamaged_vs_module_by_index.tsv"))
brown_vs_nondmg <- fread(file.path(DIRS$main, "tea_variability", "tea_brown_vs_nondamaged_comparison.tsv"))
sample_meta <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
hmm_states <- fread(file.path(DIRS$main, "hmm_states_main.tsv"))
state_fingerprints <- fread(file.path(DIRS$main, "state_fingerprints_main.tsv"))

state_key <- hmm_states[
  core == "GeoB25202_R1",
  .(mean_sst = mean(sst, na.rm = TRUE)),
  by = state
][order(mean_sst)]
state_key[, state_label := sprintf("state_%02d", seq_len(.N))]
state_key <- state_key[, .(state, state_label)]
setorder(state_key, state)

apply_state_labels <- function(dt, state_col = "state", label_col = "state_label") {
  had_label_col <- label_col %in% names(dt)
  out <- merge(dt, state_key, by.x = state_col, by.y = "state", all.x = TRUE, suffixes = c("", "_mapped"))
  if ("state_label_mapped" %in% names(out)) {
    if (had_label_col) {
      out[, state_label_source := get(label_col)]
    }
    out[, (label_col) := state_label_mapped]
    out[, state_label_mapped := NULL]
  } else if (!had_label_col) {
    out[, (label_col) := state_label]
    out[, state_label := NULL]
  }
  out
}

core_levels <- c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")
module_levels <- sort(unique(tea_sample_module$module))
module_color_map <- c(
  brown = "#8C510A",
  blue = "#1F78B4",
  green = "#33A02C",
  yellow = "#E6AB02",
  turquoise = "#1FBFCF",
  grey = "#7F7F7F",
  non_damaged = "#000000"
)

save_plot <- function(plot_obj, stem, width = 10, height = 7, dpi = 300) {
  png_file <- file.path(OUT_DIR, paste0(stem, ".png"))
  pdf_file <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  ggsave(filename = png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(filename = pdf_file, plot = plot_obj, width = width, height = height, device = cairo_pdf)
  c(png_file, pdf_file)
}

fig_index <- data.table(
  figure_id = character(),
  category = character(),
  file_png = character(),
  file_pdf = character(),
  source_table = character(),
  description = character()
)

append_index <- function(figure_id, category, saved_paths, source_table, description) {
  fig_index <<- rbind(
    fig_index,
    data.table(
      figure_id = figure_id,
      category = category,
      file_png = saved_paths[1],
      file_pdf = saved_paths[2],
      source_table = source_table,
      description = description
    ),
    fill = TRUE
  )
}

within_index_z <- function(dt, value_col, out_col = "value_scaled_by_index", index_col = "index") {
  dt_scaled <- copy(dt)
  dt_scaled[, (out_col) := {
    x <- get(value_col)
    mu <- mean(x, na.rm = TRUE)
    sig <- sd(x, na.rm = TRUE)
    if (!is.finite(sig) || sig == 0) {
      rep(0, .N)
    } else {
      (x - mu) / sig
    }
  }, by = index_col]
  dt_scaled[!is.finite(get(out_col)), (out_col) := 0]
  dt_scaled
}

# -----------------------------------------------------------------------------
# A. Module-focused TEA figures
# -----------------------------------------------------------------------------

module_long <- melt(
  tea_sample_module,
  id.vars = c("sample", "module", "module_abundance", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)[is.finite(value)]

module_long[, index := factor(index, levels = INDEX_COLS)]
module_long[, module := factor(module, levels = module_levels)]

nondmg_long <- melt(
  nondmg_sample,
  id.vars = c("sample", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)[is.finite(value)]
nondmg_long[, module := factor("non_damaged", levels = c(module_levels, "non_damaged"))]
nondmg_long[, index := factor(index, levels = INDEX_COLS)]

module_dist_long <- rbind(
  module_long,
  nondmg_long[, .(sample, module, core, age_kyr, state, state_label, index, value)],
  fill = TRUE
)
module_dist_long[, module := factor(as.character(module), levels = c(module_levels, "non_damaged"))]

module_central_raw <- module_long[, .(
  n = .N,
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE)
), by = .(module, index)]

module_central_scaled <- within_index_z(module_central_raw, "median", "median_scaled_by_index")
module_central_mean_scaled <- within_index_z(module_central_raw, "mean", "mean_scaled_by_index")

module_centered_use <- tea_centered[analysis_set == "include_grey" & is.finite(centered)]
module_centered_use[, index := factor(index, levels = INDEX_COLS)]
module_centered_use[, module := factor(module, levels = module_levels)]

module_central_centered <- module_centered_use[, .(
  n = .N,
  centered_mean = mean(centered, na.rm = TRUE),
  centered_median = median(centered, na.rm = TRUE)
), by = .(module, index)]

module_variability <- module_long[, .(
  n = .N,
  sd = sd(value, na.rm = TRUE),
  iqr = IQR(value, na.rm = TRUE)
), by = .(module, index)]

fwrite(module_central_raw, file.path(OUT_TBL, "module_central_tendency_raw.tsv"), sep = "\t")
fwrite(module_central_scaled, file.path(OUT_TBL, "module_central_tendency_raw_scaled_by_index.tsv"), sep = "\t")
fwrite(module_central_mean_scaled, file.path(OUT_TBL, "module_central_tendency_mean_scaled_by_index.tsv"), sep = "\t")
fwrite(module_central_centered, file.path(OUT_TBL, "module_central_tendency_centered_include_grey.tsv"), sep = "\t")
fwrite(module_variability, file.path(OUT_TBL, "module_variability_raw.tsv"), sep = "\t")
fwrite(module_dist_long, file.path(OUT_TBL, "module_distribution_long.tsv"), sep = "\t")

p_mod_heat_raw <- ggplot(module_central_raw, aes(x = index, y = module, fill = median)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Module-level TEA central tendency (raw)",
    subtitle = "Median TEA per module × index",
    x = "TEA index", y = "Module", fill = "Median"
  ) +
  theme_bw(base_size = 11)

p_mod_heat_centered <- ggplot(module_central_centered, aes(x = index, y = module, fill = centered_median)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Module-level TEA central tendency (within-sample centered)",
    subtitle = "Median centered TEA (include_grey analysis set)",
    x = "TEA index", y = "Module", fill = "Centered\nmedian"
  ) +
  theme_bw(base_size = 11)

p_mod_heat_raw_scaled <- ggplot(module_central_scaled, aes(x = index, y = module, fill = median_scaled_by_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Module-level TEA central tendency (raw, scaled within index)",
    subtitle = "Z-score of module medians within each index (relative enrichment/depletion)",
    x = "TEA index", y = "Module", fill = "Within-index\nz-score"
  ) +
  theme_bw(base_size = 11)

p_mod_heat_mean_scaled <- ggplot(module_central_mean_scaled, aes(x = index, y = module, fill = mean_scaled_by_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Module-level TEA central tendency (mean, scaled within index)",
    subtitle = "Z-score of module means within each index (sensitive to residual unsaturation/tails)",
    x = "TEA index", y = "Module", fill = "Within-index\nz-score"
  ) +
  theme_bw(base_size = 11)

p_mod_heat_var <- ggplot(module_variability, aes(x = index, y = module, fill = iqr)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
  labs(
    title = "Module-level TEA variability",
    subtitle = "IQR TEA per module × index",
    x = "TEA index", y = "Module", fill = "IQR"
  ) +
  theme_bw(base_size = 11)

p_mod_dist <- ggplot(module_dist_long, aes(x = module, y = value, fill = module)) +
  geom_boxplot(outlier.alpha = 0.15, width = 0.75) +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = module_color_map, drop = FALSE) +
  guides(fill = "none") +
  labs(
    title = "Distribution of module-level TEA values",
    subtitle = "Sample × module values across indices",
    x = "Module", y = "TEA value"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

append_index("A1", "module", save_plot(p_mod_heat_raw, "tea_module_heatmap_central_raw", 8, 5.5),
             "tables/module_central_tendency_raw.tsv",
             "Heatmap of module-level raw TEA medians across indices.")
append_index("A1c", "module", save_plot(p_mod_heat_raw_scaled, "tea_module_heatmap_central_raw_scaled_by_index", 8, 5.5),
             "tables/module_central_tendency_raw_scaled_by_index.tsv",
             "Within-index z-scored module TEA medians for cross-index visual comparability (relative, not absolute magnitude).")
append_index("A1d", "module", save_plot(p_mod_heat_mean_scaled, "tea_module_heatmap_central_mean_scaled_by_index", 8, 5.5),
             "tables/module_central_tendency_mean_scaled_by_index.tsv",
             "Within-index z-scored module TEA means; exploratory complement to median-scaled views, useful for saturated indices with median collapse.")
append_index("A1b", "module", save_plot(p_mod_heat_centered, "tea_module_heatmap_central_centered_include_grey", 8, 5.5),
             "tables/module_central_tendency_centered_include_grey.tsv",
             "Heatmap of within-sample centered module TEA medians (include_grey set).")
append_index("A2", "module", save_plot(p_mod_heat_var, "tea_module_heatmap_variability_iqr", 8, 5.5),
             "tables/module_variability_raw.tsv",
             "Heatmap of module-level TEA IQR across indices.")
append_index("A3", "module", save_plot(p_mod_dist, "tea_module_distribution_boxplot_by_index", 11, 7.5),
             "tables/module_distribution_long.tsv",
             "Faceted module TEA distributions by index.")

# -----------------------------------------------------------------------------
# B. State-focused TEA figures
# -----------------------------------------------------------------------------

sample_long <- melt(
  tea_sample,
  id.vars = c("sample", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)[is.finite(value)]

sample_long <- apply_state_labels(sample_long)

sample_long[, index := factor(index, levels = INDEX_COLS)]
sample_long[, state_label := factor(state_label, levels = paste0("state_", sprintf("%02d", seq_len(nrow(state_key)))))]

state_central <- sample_long[, .(
  n = .N,
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE)
), by = .(state, state_label, index)][order(state, index)]

state_central_scaled <- within_index_z(state_central, "median", "median_scaled_by_index")
state_central_mean_scaled <- within_index_z(state_central, "mean", "mean_scaled_by_index")

state_module_plot <- copy(state_module_specialization)
state_module_plot <- state_module_plot[is.finite(centered_median)]
state_module_plot <- apply_state_labels(state_module_plot)
state_module_plot[, index := factor(index, levels = INDEX_COLS)]
state_module_plot[, module := factor(module, levels = module_levels)]
state_module_plot[, state_label := factor(state_label, levels = paste0("state_", sprintf("%02d", seq_len(nrow(state_key)))))]

fwrite(state_central, file.path(OUT_TBL, "state_central_tendency_raw.tsv"), sep = "\t")
fwrite(state_central_scaled, file.path(OUT_TBL, "state_central_tendency_raw_scaled_by_index.tsv"), sep = "\t")
fwrite(state_central_mean_scaled, file.path(OUT_TBL, "state_central_tendency_mean_scaled_by_index.tsv"), sep = "\t")
fwrite(state_module_plot, file.path(OUT_TBL, "state_module_specialization_centered_median.tsv"), sep = "\t")
fwrite(sample_long, file.path(OUT_TBL, "state_distribution_long.tsv"), sep = "\t")

p_state_heat <- ggplot(state_central, aes(x = index, y = state_label, fill = median)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "State-level TEA central tendency",
    subtitle = "Median sample-level TEA per state × index",
    x = "TEA index", y = "State", fill = "Median"
  ) +
  theme_bw(base_size = 11)

p_state_mod_spec <- ggplot(state_module_plot, aes(x = module, y = state_label, fill = centered_median)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~index, ncol = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "State × module specialization (centered TEA)",
    subtitle = "Centered median values from 04b specialization summary",
    x = "Module", y = "State", fill = "Centered\nmedian"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_state_heat_scaled <- ggplot(state_central_scaled, aes(x = index, y = state_label, fill = median_scaled_by_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "State-level TEA central tendency (raw, scaled within index)",
    subtitle = "Z-score of state medians within each index (relative enrichment/depletion)",
    x = "TEA index", y = "State", fill = "Within-index\nz-score"
  ) +
  theme_bw(base_size = 11)

p_state_heat_mean_scaled <- ggplot(state_central_mean_scaled, aes(x = index, y = state_label, fill = mean_scaled_by_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "State-level TEA central tendency (mean, scaled within index)",
    subtitle = "Z-score of state means within each index (sensitive to residual unsaturation/tails)",
    x = "TEA index", y = "State", fill = "Within-index\nz-score"
  ) +
  theme_bw(base_size = 11)

p_state_dist <- ggplot(sample_long, aes(x = state_label, y = value, fill = state_label)) +
  geom_boxplot(outlier.alpha = 0.15, width = 0.75) +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  guides(fill = "none") +
  labs(
    title = "Sample-level TEA distributions by HMM state",
    x = "State", y = "TEA value"
  ) +
  theme_bw(base_size = 10)

append_index("B4", "state", save_plot(p_state_heat, "tea_state_heatmap_central_raw", 8, 5),
             "tables/state_central_tendency_raw.tsv",
             "Heatmap of state-level sample TEA medians across indices.")
append_index("B4b", "state", save_plot(p_state_heat_scaled, "tea_state_heatmap_central_raw_scaled_by_index", 8, 5),
             "tables/state_central_tendency_raw_scaled_by_index.tsv",
             "Within-index z-scored state TEA medians for cross-index visual comparability (relative, not absolute magnitude).")
append_index("B4c", "state", save_plot(p_state_heat_mean_scaled, "tea_state_heatmap_central_mean_scaled_by_index", 8, 5),
             "tables/state_central_tendency_mean_scaled_by_index.tsv",
             "Within-index z-scored state TEA means; exploratory complement to median-scaled views, useful for saturated indices with median collapse.")
append_index("B5", "state", save_plot(p_state_mod_spec, "tea_state_module_specialization_heatmap", 12, 8),
             "tables/state_module_specialization_centered_median.tsv",
             "State × module centered-median specialization heatmap across indices.")
append_index("B6", "state", save_plot(p_state_dist, "tea_state_distribution_boxplot_by_index", 10.5, 7.5),
             "tables/state_distribution_long.tsv",
             "Sample-level TEA distributions grouped by HMM state.")

# State fingerprint (eigengene profile) heatmaps
state_fingerprint_plot <- apply_state_labels(copy(state_fingerprints))
me_cols <- grep("^ME", names(state_fingerprint_plot), value = TRUE)
if (length(me_cols) == 0) {
  stop("state_fingerprints_main.tsv has no module eigengene columns (expected names starting with 'ME').")
}

state_fingerprint_long <- melt(
  state_fingerprint_plot,
  id.vars = intersect(c("state", "state_label", "n_samples", "mean_age_kyr"), names(state_fingerprint_plot)),
  measure.vars = me_cols,
  variable.name = "module",
  value.name = "eigengene_mean"
)[is.finite(eigengene_mean)]

state_fingerprint_long[, state_label := factor(state_label, levels = paste0("state_", sprintf("%02d", seq_len(nrow(state_key)))))]
state_fingerprint_long[, module := factor(module, levels = me_cols)]

state_fingerprint_long_scaled <- within_index_z(
  state_fingerprint_long,
  value_col = "eigengene_mean",
  out_col = "eigengene_z_by_module",
  index_col = "module"
)

fwrite(state_fingerprint_long, file.path(OUT_TBL, "state_fingerprint_module_eigengene_raw.tsv"), sep = "\t")
fwrite(state_fingerprint_long_scaled, file.path(OUT_TBL, "state_fingerprint_module_eigengene_scaled_by_module.tsv"), sep = "\t")

p_state_fingerprint_scaled <- ggplot(state_fingerprint_long_scaled, aes(x = module, y = state_label, fill = eigengene_z_by_module)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "HMM state fingerprints across module eigengenes",
    subtitle = "Within-module z-score across states (relative high/low loading per module)",
    x = "Module eigengene", y = "State", fill = "Within-module\nz-score"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_state_fingerprint_raw <- ggplot(state_fingerprint_long, aes(x = module, y = state_label, fill = eigengene_mean)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "HMM state fingerprints across module eigengenes (raw)",
    subtitle = "Mean eigengene value per state",
    x = "Module eigengene", y = "State", fill = "Mean\neigengene"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

append_index("B7", "state", save_plot(p_state_fingerprint_scaled, "state_fingerprints_module_heatmap", 9, 4.8),
             "tables/state_fingerprint_module_eigengene_scaled_by_module.tsv",
             "State × module eigengene fingerprint heatmap scaled within module; shows relative high/low loading across states, not absolute across-module magnitude.")
append_index("B7b", "state", save_plot(p_state_fingerprint_raw, "state_fingerprints_module_heatmap_raw", 9, 4.8),
             "tables/state_fingerprint_module_eigengene_raw.tsv",
             "State × module eigengene fingerprint heatmap using raw mean eigengene values per state.")

# -----------------------------------------------------------------------------
# C. Time-focused TEA figures
# -----------------------------------------------------------------------------

sample_time <- copy(sample_long)
sample_time[, core := factor(core, levels = core_levels)]

module_plus_nondmg_time <- rbind(
  module_long[, .(sample, module, core, age_kyr, state, state_label, index, value)],
  melt(
    nondmg_sample,
    id.vars = c("sample", "core", "age_kyr", "state", "state_label"),
    measure.vars = INDEX_COLS,
    variable.name = "index",
    value.name = "value"
  )[is.finite(value), .(sample, module = "non_damaged", core, age_kyr, state, state_label, index, value)],
  fill = TRUE
)
module_plus_nondmg_time <- apply_state_labels(module_plus_nondmg_time)
module_plus_nondmg_time[, core := factor(core, levels = core_levels)]
module_plus_nondmg_time[, index := factor(index, levels = INDEX_COLS)]
module_plus_nondmg_time[, module := factor(module, levels = c("non_damaged", module_levels))]

module_time <- merge(
  module_centered_use[, .(sample, module, core, age_kyr, state, state_label, index, centered)],
  state_key,
  by = "state",
  all.x = TRUE,
  suffixes = c("", "_ref")
)
module_time <- apply_state_labels(module_time)
module_time[, core := factor(core, levels = core_levels)]
module_time[, module := factor(module, levels = module_levels)]
module_time[, index := factor(index, levels = INDEX_COLS)]

fwrite(sample_time, file.path(OUT_TBL, "sample_time_long.tsv"), sep = "\t")
fwrite(module_time, file.path(OUT_TBL, "module_time_centered_long.tsv"), sep = "\t")
fwrite(module_plus_nondmg_time, file.path(OUT_TBL, "module_plus_nondamaged_time_long.tsv"), sep = "\t")
fwrite(temporal_summary, file.path(OUT_TBL, "temporal_descriptive_summary_input_copy.tsv"), sep = "\t")

p_time_sample <- ggplot(sample_time, aes(x = age_kyr, y = value, color = core)) +
  geom_point(alpha = 0.45, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.6, span = 0.6) +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  scale_x_reverse(labels = label_number(accuracy = 1)) +
  labs(
    title = "Sample-level TEA through time across cores",
    subtitle = "Exploratory scatter + loess smooth; colors denote core",
    x = "Age (kyr)", y = "TEA value", color = "Core"
  ) +
  theme_bw(base_size = 10)

p_time_module <- ggplot(module_time, aes(x = age_kyr, y = centered, color = core)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey35") +
  geom_point(alpha = 0.35, size = 0.75) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.45, span = 0.6) +
  facet_grid(index ~ module, scales = "free_y") +
  scale_x_reverse(labels = label_number(accuracy = 1)) +
  labs(
    title = "Module-level centered TEA through time",
    subtitle = "Within-sample centered TEA (include_grey), faceted by index × module",
    x = "Age (kyr)", y = "Centered TEA", color = "Core"
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")

p_time_module_plus_nondmg <- ggplot(module_plus_nondmg_time, aes(x = age_kyr, y = value, color = module)) +
  geom_point(alpha = 0.35, size = 0.75) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.45, span = 0.7) +
  facet_grid(index ~ core, scales = "free_y") +
  scale_color_manual(values = module_color_map) +
  scale_x_reverse(labels = label_number(accuracy = 1)) +
  labs(
    title = "Module + non-damaged TEA through time",
    subtitle = "Raw sample-level TEA values; non_damaged is a pseudo-module for comparison",
    x = "Age (kyr)", y = "TEA value", color = "Module /\ncomparison"
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")

append_index("C7", "time", save_plot(p_time_sample, "tea_sample_time_scatter_smooth_by_index", 11, 7.5),
             "tables/sample_time_long.tsv",
             "Sample-level TEA through time across all cores, faceted by index.")
append_index("C8", "time", save_plot(p_time_module, "tea_module_centered_time_scatter_smooth_grid", 16, 10),
             "tables/module_time_centered_long.tsv",
             "Module-level centered TEA through time (index × module grid).")
append_index("C9", "time", save_plot(p_time_module_plus_nondmg, "tea_module_plus_nondamaged_time_scatter_smooth_grid", 12, 10),
             "tables/module_plus_nondamaged_time_long.tsv",
             "Raw TEA through time by core × index for modules plus pseudo-module non_damaged.")

# -----------------------------------------------------------------------------
# D. Optional comparison figures (non-damaged focus)
# -----------------------------------------------------------------------------

nondmg_profile <- copy(nondmg_vs_module_by_index)
nondmg_profile[, index := factor(index, levels = INDEX_COLS)]
nondmg_profile[, module := factor(module, levels = module_levels)]

nondmg_pseudomod <- unique(nondmg_profile[, .(
  module = "non_damaged",
  index,
  profile_mean = nondmg_mean
)])

module_profile <- nondmg_profile[, .(module, index, profile_mean = module_mean)]
module_plus_nd <- rbind(module_profile, nondmg_pseudomod, fill = TRUE)
module_plus_nd[, module := factor(module, levels = c("non_damaged", module_levels))]
module_plus_nd_scaled <- within_index_z(module_plus_nd, "profile_mean", "profile_mean_scaled_by_index")

brown_focus <- brown_vs_nondmg[comparison_scope == "per_index_aggregate" & is.finite(diff_mean_brown_minus_nondmg)]
brown_focus[, index := factor(index, levels = INDEX_COLS)]

fwrite(module_plus_nd, file.path(OUT_TBL, "nondamaged_plus_module_profile_mean.tsv"), sep = "\t")
fwrite(module_plus_nd_scaled, file.path(OUT_TBL, "nondamaged_plus_module_profile_mean_scaled_by_index.tsv"), sep = "\t")
fwrite(brown_focus, file.path(OUT_TBL, "brown_vs_nondamaged_per_index.tsv"), sep = "\t")
fwrite(dominance_calls, file.path(OUT_TBL, "dominance_calls_input_copy.tsv"), sep = "\t")
fwrite(nondmg_vs_module, file.path(OUT_TBL, "nondamaged_vs_module_comparison_input_copy.tsv"), sep = "\t")

p_nd_heat <- ggplot(module_plus_nd, aes(x = index, y = module, fill = profile_mean)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Non-damaged aggregate TEA profile vs damaged modules",
    subtitle = "Pseudo-module 'non_damaged' included for direct profile comparison",
    x = "TEA index", y = "Profile row", fill = "Mean"
  ) +
  theme_bw(base_size = 11)

p_brown_nd <- ggplot(brown_focus, aes(x = index)) +
  geom_segment(aes(xend = index, y = nondmg_mean, yend = brown_mean), color = "grey40", linewidth = 1) +
  geom_point(aes(y = nondmg_mean, color = "Non-damaged"), size = 2.2) +
  geom_point(aes(y = brown_mean, color = "Brown module"), size = 2.2) +
  scale_color_manual(values = c("Non-damaged" = "#2166ac", "Brown module" = "#b2182b")) +
  labs(
    title = "Brown module vs non-damaged aggregate TEA",
    subtitle = "Per-index mean comparison (segment links values)",
    x = "TEA index", y = "Mean TEA", color = NULL
  ) +
  theme_bw(base_size = 11)

p_nd_heat_scaled <- ggplot(module_plus_nd_scaled, aes(x = index, y = module, fill = profile_mean_scaled_by_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Non-damaged aggregate TEA profile vs damaged modules (scaled within index)",
    subtitle = "Z-score across rows within each index (relative enrichment/depletion)",
    x = "TEA index", y = "Profile row", fill = "Within-index\nz-score"
  ) +
  theme_bw(base_size = 11)

p_nd_heat_mean_scaled <- ggplot(module_plus_nd_scaled, aes(x = index, y = module, fill = profile_mean_scaled_by_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Non-damaged aggregate TEA profile vs damaged modules (mean, scaled within index)",
    subtitle = "Within-index z-score of profile means (sensitive to residual unsaturation/tails)",
    x = "TEA index", y = "Profile row", fill = "Within-index\nz-score"
  ) +
  theme_bw(base_size = 11)

append_index("D9", "comparison", save_plot(p_nd_heat, "tea_nondamaged_vs_modules_heatmap", 8.5, 5.8),
             "tables/nondamaged_plus_module_profile_mean.tsv",
             "Heatmap comparing non-damaged aggregate TEA profile to all module profiles.")
append_index("D9b", "comparison", save_plot(p_nd_heat_scaled, "tea_nondamaged_vs_modules_heatmap_scaled_by_index", 8.5, 5.8),
             "tables/nondamaged_plus_module_profile_mean_scaled_by_index.tsv",
             "Within-index z-scored non-damaged vs module profile heatmap for cross-index visual comparability (relative, not absolute magnitude).")
append_index("D9c", "comparison", save_plot(p_nd_heat_mean_scaled, "tea_nondamaged_vs_modules_heatmap_mean_scaled_by_index", 8.5, 5.8),
             "tables/nondamaged_plus_module_profile_mean_scaled_by_index.tsv",
             "Within-index z-scored mean-profile heatmap as exploratory complement where saturated indices retain tail-driven differences.")
append_index("D10", "comparison", save_plot(p_brown_nd, "tea_brown_vs_nondamaged_dotsegment", 8.5, 5.5),
             "tables/brown_vs_nondamaged_per_index.tsv",
             "Per-index dot/segment comparison between brown module and non-damaged profile.")

# -----------------------------------------------------------------------------
# Output hygiene
# -----------------------------------------------------------------------------

fig_index <- fig_index[order(figure_id)]
fwrite(fig_index, file.path(OUT_DIR, "tea_figures_index.tsv"), sep = "\t")

index_md <- c(
  "# TEA exploratory figure index",
  "",
  "Generated by `code/wgcna_hmm/08_tea_figures.R`.",
  "",
  "| Figure ID | Category | PNG | PDF | Source table | Description |",
  "|---|---|---|---|---|---|"
)

index_rows <- fig_index[, sprintf(
  "| %s | %s | %s | %s | %s | %s |",
  figure_id, category,
  basename(file_png), basename(file_pdf), source_table, description
)]

writeLines(c(index_md, index_rows), con = file.path(OUT_DIR, "tea_figures_index.md"))

notes <- data.table(
  note_key = c(
    "analysis_scope",
    "index_set",
    "module_time_view",
    "module_plus_nondamaged_time_view",
    "state_labels",
    "interpretation_guardrail",
    "scaled_heatmap_interpretation_guardrail",
    "median_vs_mean_scaled_heatmaps",
    "mean_scaled_heatmap_interpretation_guardrail",
    "assumption_age_axis",
    "state_fingerprint_interpretation_guardrail"
  ),
  note_value = c(
    "Exploratory TEA figure battery using already-generated clean main-window outputs.",
    paste(INDEX_COLS, collapse = ", "),
    "Module-through-time plot uses within-sample centered TEA (analysis_set = include_grey) to emphasize specialization.",
    "Additional module/non_damaged time plot uses raw TEA values; non_damaged is a pseudo-module row used only for visual comparison.",
    "State labels follow 07_figures.R logic: states ordered by GeoB25202_R1 mean SST and relabeled state_01..state_05.",
    "TEA reflects annotation-weighted metabolic potential, not direct measured activity.",
    "Scaled-by-index heatmaps use within-index z-scores computed across displayed rows (modules/states/profile rows). They indicate relative enrichment/depletion within each index and are not comparable as absolute TEA magnitude across indices.",
    "Median-scaled heatmaps emphasize robust central tendency and resist outliers; mean-scaled heatmaps are added as complementary exploratory views.",
    "Mean-scaled heatmaps are more sensitive to small tail differences and partial unsaturation (e.g., MGI/MFI) and should be interpreted as supportive/exploratory rather than standalone evidence.",
    "All time panels use age_kyr as provided in workflow outputs; x-axis reversed for paleotime readability.",
    "State fingerprint heatmaps summarize module eigengene means by HMM state (not TEA). The scaled version z-scores each module across states, so color indicates relative high/low state loading within module, not absolute eigengene magnitude across modules."
  )
)

fwrite(notes, file.path(OUT_DIR, "tea_figures_notes.tsv"), sep = "\t")

write_run_metadata(
  file.path(OUT_DIR, "08_tea_figures_run_metadata.tsv"),
  "08_tea_figures.R",
  extra = list(
    n_figures = nrow(fig_index),
    output_dir = OUT_DIR,
    n_samples = uniqueN(tea_sample$sample),
    n_sample_module_rows = nrow(tea_sample_module),
    n_centered_rows = nrow(module_centered_use),
    interpretation_scope = "exploratory_tea_potential_not_activity"
  )
)
write_session_info(file.path(OUT_DIR, "08_tea_figures_sessionInfo.txt"))

log_msg("08_tea_figures complete.")
