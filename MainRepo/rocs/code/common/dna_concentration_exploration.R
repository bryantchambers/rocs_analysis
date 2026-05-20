#!/usr/bin/env Rscript

# Curated supplementary pipeline for library_concentration vs proxy analyses.
#
# Goal:
# - Recompute core analyses from manuscript-vetted stage scripts.
# - Integrate only the outputs needed for Supplementary Results narrative.
# - Write a concise final package under:
#     results/common/proxies/final/
#     plots/proxies/final/
#
# Guardrails:
# - observational associations only (not causal proof)
# - no downstream sequencing-output metrics as explanatory predictors
# - no broad directory copy of all stage outputs

options(stringsAsFactors = FALSE)

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

safe_read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required file missing: ", path)
  utils::read.delim(path, check.names = FALSE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Required file missing: ", path)
  utils::read.csv(path, check.names = FALSE)
}

safe_spearman <- function(x, y, min_n = 12) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < min_n || length(unique(x2)) < 3 || length(unique(y2)) < 3) {
    return(list(n = length(x2), rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(n = length(x2), rho = unname(ct$estimate), p = ct$p.value)
}

nearest_join_by_depth <- function(reference_df, reference_depth_col, source_df, source_depth_col, tol_cm = 3) {
  ref_depth <- as.numeric(reference_df[[reference_depth_col]])
  src_depth <- as.numeric(source_df[[source_depth_col]])

  out <- reference_df
  src_cols <- setdiff(names(source_df), source_depth_col)
  for (nm in src_cols) out[[nm]] <- NA
  out$proxy_matched_depth_cm <- NA_real_
  out$proxy_depth_diff_cm <- NA_real_
  out$proxy_match_within_tolerance <- FALSE

  for (i in seq_along(ref_depth)) {
    if (!is.finite(ref_depth[i])) next
    j <- which.min(abs(src_depth - ref_depth[i]))
    if (length(j) == 0 || !is.finite(src_depth[j])) next
    dd <- abs(src_depth[j] - ref_depth[i])
    out$proxy_matched_depth_cm[i] <- src_depth[j]
    out$proxy_depth_diff_cm[i] <- dd
    if (is.finite(dd) && dd <= tol_cm) {
      out$proxy_match_within_tolerance[i] <- TRUE
      for (nm in src_cols) out[[nm]][i] <- source_df[[nm]][j]
    }
  }
  out
}

run_r_script <- function(script_path, log_file) {
  result <- tryCatch(
    {
      x <- system2("Rscript", args = script_path, stdout = TRUE, stderr = TRUE)
      st <- attr(x, "status")
      if (is.null(st)) st <- 0L
      list(stdout = x, status = as.integer(st), error = NA_character_)
    },
    error = function(e) {
      list(stdout = paste("system2 error:", conditionMessage(e)), status = 999L, error = conditionMessage(e))
    }
  )

  log_lines <- c(
    paste0("timestamp_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("script\t", script_path),
    paste0("status\t", result$status),
    "---- begin stdout_stderr ----",
    result$stdout,
    "---- end stdout_stderr ----"
  )
  writeLines(log_lines, con = log_file)

  result
}

write_geob_proxy_plot <- function(df, out_png, out_pdf) {
  if (nrow(df) == 0) return(invisible(NULL))
  ord <- order(abs(df$rho_age_adjusted_linear), decreasing = TRUE, na.last = NA)
  plot_df <- df[ord, , drop = FALSE]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  draw_plot <- function(device_fun, path) {
    device_fun(path)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mar = c(8, 5, 3, 1))
    vals <- plot_df$rho_age_adjusted_linear
    cols <- ifelse(vals < 0, "#2C7BB6", "#D7191C")
    mids <- barplot(
      height = vals,
      names.arg = plot_df$proxy,
      las = 2,
      col = cols,
      border = NA,
      ylab = "Age-adjusted Spearman rho",
      main = "GeoB25202 extra proxies vs library concentration"
    )
    abline(h = 0, lty = 2, col = "grey40")
    text(
      x = mids,
      y = ifelse(vals >= 0, vals + 0.02, vals - 0.02),
      labels = sprintf("n=%d", plot_df$n),
      cex = 0.75,
      pos = ifelse(vals >= 0, 3, 1)
    )
    dev.off()
  }

  draw_plot(function(p) grDevices::png(p, width = 1400, height = 900, res = 150), out_png)
  draw_plot(function(p) grDevices::pdf(p, width = 11, height = 7), out_pdf)
}

repo_required_files <- c(
  file.path("data", "metadata_v5.tsv"),
  file.path("data", "combined_xrf_geochemistry_curated.csv"),
  file.path("data", "combined_foraminifera_geochem.tsv"),
  file.path("data", "combined_sst_proxies_separate_columns.tsv"),
  file.path("data", "geob25202_clean_proxies.tsv")
)

missing_required <- repo_required_files[!file.exists(repo_required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input files (run from repository root): ", paste(missing_required, collapse = ", "))
}

results_root <- file.path("results", "common", "proxies")
plots_root <- file.path("plots", "proxies")
final_results <- file.path(results_root, "final")
final_plots <- file.path(plots_root, "final")
logs_dir <- file.path(results_root, "logs")

dir_create(results_root)
dir_create(plots_root)
dir_create(final_results)
dir_create(final_plots)
dir_create(logs_dir)

for (f in list.files(final_results, full.names = TRUE, recursive = TRUE, include.dirs = FALSE)) file.remove(f)
for (f in list.files(final_plots, full.names = TRUE, recursive = TRUE, include.dirs = FALSE)) file.remove(f)

stage_df <- data.frame(
  stage_id = c("exploratory", "st13", "grouped_axes", "state_dependent", "multicore"),
  script_path = c(
    file.path("code", "common", "explore_library_concentration_vs_proxies.R"),
    file.path("code", "common", "analyze_st13_library_concentration_model.R"),
    file.path("code", "common", "analyze_st13_proxy_axes_model.R"),
    file.path("code", "common", "analyze_st13_state_dependent_terrigenous_model.R"),
    file.path("code", "common", "analyze_multicore_rbsr_validation.R")
  ),
  stringsAsFactors = FALSE
)

if (!all(file.exists(stage_df$script_path))) {
  stop("Missing stage script(s): ", paste(stage_df$script_path[!file.exists(stage_df$script_path)], collapse = ", "))
}

run_manifest <- data.frame(
  stage_id = stage_df$stage_id,
  script_path = stage_df$script_path,
  status = NA_character_,
  run_status_code = NA_integer_,
  log_file = file.path(logs_dir, paste0(stage_df$stage_id, ".log")),
  stringsAsFactors = FALSE
)

message("Running prerequisite analyses (no wholesale output copy)...")
for (i in seq_len(nrow(stage_df))) {
  sid <- stage_df$stage_id[i]
  spath <- stage_df$script_path[i]
  lpath <- run_manifest$log_file[i]
  message(sprintf("[%d/%d] %s", i, nrow(stage_df), spath))
  run <- run_r_script(spath, lpath)
  run_manifest$run_status_code[i] <- as.integer(run$status)
  if (!identical(as.integer(run$status), 0L)) {
    run_manifest$status[i] <- "failed"
    write_tsv(run_manifest, file.path(final_results, "pipeline_run_manifest.tsv"))
    stop("Stage failed: ", sid, " (status=", run$status, "). See log: ", lpath)
  }
  run_manifest$status[i] <- "success"
}

write_tsv(run_manifest, file.path(final_results, "pipeline_run_manifest.tsv"))

# ---- Load stage outputs needed for curated integration ----
exploratory_overall <- safe_read_tsv(file.path("analysis", "library_concentration_proxy_exploration", "correlations_overall.tsv"))

st13_layer <- safe_read_tsv(file.path("analysis", "st13_library_concentration_model", "st13_layer_level_analysis_table.tsv"))
st13_screen <- safe_read_tsv(file.path("analysis", "st13_library_concentration_model", "screening_correlations_raw_and_detrended.tsv"))
st13_primary <- safe_read_tsv(file.path("analysis", "st13_library_concentration_model", "primary_model_summary_table.tsv"))
st13_auto <- safe_read_tsv(file.path("analysis", "st13_library_concentration_model", "autocorrelation_robust_model_summary.tsv"))

axis_membership <- safe_read_tsv(file.path("analysis", "st13_proxy_axes_model", "axis_variable_membership.tsv"))
axis_loadings <- safe_read_tsv(file.path("analysis", "st13_proxy_axes_model", "axis_loadings_or_weights.tsv"))
axis_grouped <- safe_read_tsv(file.path("analysis", "st13_proxy_axes_model", "grouped_axis_model_summary.tsv"))
axis_nested <- safe_read_tsv(file.path("analysis", "st13_proxy_axes_model", "nested_axis_model_comparison.tsv"))

state_interaction <- safe_read_tsv(file.path("analysis", "st13_state_dependent_terrigenous_model", "interaction_model_summary.tsv"))
state_comparison <- safe_read_tsv(file.path("analysis", "st13_state_dependent_terrigenous_model", "interaction_model_comparison.tsv"))

multi_summary <- safe_read_tsv(file.path("analysis", "multicore_rbsr_library_validation", "multicore_rbsr_summary.tsv"))
multi_detail <- safe_read_tsv(file.path("analysis", "multicore_rbsr_library_validation", "multicore_rbsr_effect_details.tsv"))

# ---- A. Exploratory candidate summary (concise) ----
preselected_candidates <- c(
  "rb_sr", "rb", "zr", "al_si", "ti_al", "si_ti", "ca_ti", "ba_ti", "mn_fe",
  "sum_greenland", "sum_illites_micas", "sum_montmorillonites_smectites", "sum_expandable_clays",
  "toc_percent_weight", "ratio_greenland_to_iceland"
)

exploratory_candidates <- exploratory_overall[
  exploratory_overall$n_matched >= 50 &
    (
      exploratory_overall$proxy_name %in% preselected_candidates |
        (is.finite(exploratory_overall$q_value_bh) & exploratory_overall$q_value_bh <= 0.02 & abs(exploratory_overall$spearman_rho) >= 0.30)
    ),
  c("source_dataset", "proxy_name", "n_matched", "spearman_rho", "spearman_p", "q_value_bh", "library_vs_age_spearman", "proxy_vs_age_spearman")
]

exploratory_candidates <- exploratory_candidates[order(exploratory_candidates$q_value_bh, -abs(exploratory_candidates$spearman_rho)), , drop = FALSE]
if (nrow(exploratory_candidates) > 18) exploratory_candidates <- exploratory_candidates[seq_len(18), , drop = FALSE]

exploratory_candidates$candidate_group <- ifelse(
  exploratory_candidates$proxy_name %in% c("rb_sr", "rb", "zr", "al_si", "ti_al", "si_ti"), "terrigenous/mineralogical",
  ifelse(
    exploratory_candidates$proxy_name %in% c("sum_greenland", "sum_illites_micas", "sum_montmorillonites_smectites", "sum_expandable_clays"),
    "GeoB25202 clay/mineral",
    ifelse(exploratory_candidates$proxy_name %in% c("toc_percent_weight"), "GeoB25202 TOC", "other")
  )
)

write_tsv(exploratory_candidates, file.path(final_results, "exploratory_candidate_proxy_summary.tsv"))

# ---- B. ST13 primary evidence ----
st13_essential_cols <- c(
  "depth_in_core_cm", "layer_id", "y_bp", "library_concentration", "library_concentration_log10",
  "rb_sr", "rb", "zr", "ba_ti", "ca_ti", "mn_fe", "n_xrf_rows_median"
)
st13_essential_cols <- st13_essential_cols[st13_essential_cols %in% names(st13_layer)]
st13_essential <- st13_layer[, st13_essential_cols, drop = FALSE]
write_tsv(st13_essential, file.path(final_results, "st13_layer_table_essential.tsv"))

st13_screen_keep <- c(
  "rb_sr", "rb", "zr", "al_si", "ti_al", "si_ti",
  "ba_ti", "ca_ti", "ca", "mn_fe", "fe_al", "mn", "fe",
  "foram_n_pachyderma_d13c_vpdb_corr", "foram_n_pachyderma_d18o_vpdb_corr"
)
st13_screen_curated <- st13_screen[
  st13_screen$variable %in% st13_screen_keep & st13_screen$n_detrended >= 20,
  c("variable", "family", "n_raw", "rho_raw", "p_raw", "n_detrended", "rho_detrended", "p_detrended", "q_detrended_bh_family")
]
st13_screen_curated <- st13_screen_curated[order(abs(st13_screen_curated$rho_detrended), decreasing = TRUE), , drop = FALSE]
write_tsv(st13_screen_curated, file.path(final_results, "st13_screening_selected_proxies.tsv"))

st13_primary_curated <- st13_primary[st13_primary$term %in% c("z_rb_sr", "z_ba_ti", "z_ca_ti", "z_mn_fe"), , drop = FALSE]
st13_primary_curated$model_id <- "lm_primary"
st13_boot_curated <- st13_auto[
  st13_auto$model_id %in% c("lm_block_bootstrap_primary") &
    st13_auto$term %in% c("z_rb_sr", "z_ba_ti", "z_ca_ti", "z_mn_fe"),
  c("model_id", "term", "estimate", "std_error", "p_value", "ci_low", "ci_high")
]
st13_primary_boot <- rbind(
  st13_primary_curated[, c("model_id", "term", "estimate", "std_error", "p_value", "ci_low", "ci_high")],
  st13_boot_curated
)
write_tsv(st13_primary_boot, file.path(final_results, "st13_primary_model_and_bootstrap_summary.tsv"))

# ---- C. ST13 grouped-axis collinearity resolution ----
axis_membership_loadings <- merge(
  axis_membership[, c("axis", "variable", "included", "exclusion_reason"), drop = FALSE],
  axis_loadings[, c("axis", "variable", "oriented_loading_or_weight"), drop = FALSE],
  by = c("axis", "variable"),
  all.x = TRUE,
  all.y = FALSE
)
write_tsv(axis_membership_loadings, file.path(final_results, "st13_grouped_axis_membership_loadings.tsv"))

write_tsv(axis_grouped, file.path(final_results, "st13_grouped_axis_model_summary.tsv"))
write_tsv(axis_nested, file.path(final_results, "st13_grouped_axis_nested_comparison.tsv"))

# ---- D. ST13 state dependence (concise, negative/inconclusive) ----
state_interaction_key <- state_interaction[, c("model_id", "modifier_type", "term", "estimate", "p_value", "ci_low", "ci_high", "interaction_supported_p_lt_0_05"), drop = FALSE]
state_join <- merge(
  state_interaction_key,
  state_comparison[, c("modifier_type", "delta_aic_with_minus_without", "p_nested"), drop = FALSE],
  by = "modifier_type",
  all.x = TRUE,
  all.y = FALSE
)
write_tsv(state_join, file.path(final_results, "st13_state_dependence_interaction_summary.tsv"))

# ---- E. Multicore Rb/Sr validation ----
multi_curated <- merge(
  multi_detail[, c(
    "core", "n_obs", "rb_sr_raw_rho", "rb_sr_raw_p", "rb_sr_detrended_rho", "rb_sr_detrended_p",
    "rb_sr_model_estimate", "rb_sr_model_p", "rb_sr_model_ci_low", "rb_sr_model_ci_high", "rb_sr_model_type"
  ), drop = FALSE],
  multi_summary[, c("core", "rb_sr_supported_yes_no", "notes"), drop = FALSE],
  by = "core",
  all.x = TRUE,
  all.y = FALSE
)
multi_curated <- multi_curated[order(match(multi_curated$core, c("ST13", "ST5", "ST8", "GeoB25202"))), , drop = FALSE]
write_tsv(multi_curated, file.path(final_results, "multicore_rbsr_validation_summary.tsv"))

# ---- F. GeoB25202 extra proxies (raw + simple detrended) ----
metadata <- safe_read_tsv(file.path("data", "metadata_v5.tsv"))
geob_clean <- safe_read_tsv(file.path("data", "geob25202_clean_proxies.tsv"))

meta_geob <- metadata[grepl("^GeoB25202", metadata$core), c("depth_in_core_cm", "y_bp", "library_concentration"), drop = FALSE]
meta_geob <- meta_geob[is.finite(meta_geob$depth_in_core_cm) & is.finite(meta_geob$library_concentration), , drop = FALSE]
meta_geob$depth_in_core_cm <- as.numeric(meta_geob$depth_in_core_cm)
meta_geob$y_bp <- as.numeric(meta_geob$y_bp)
meta_geob$library_concentration <- as.numeric(meta_geob$library_concentration)

meta_geob_layer <- stats::aggregate(
  cbind(y_bp, library_concentration) ~ depth_in_core_cm,
  data = meta_geob,
  FUN = function(x) mean(x, na.rm = TRUE)
)
meta_geob_layer$library_concentration_log10 <- log10(meta_geob_layer$library_concentration)

geob_proxy <- geob_clean[geob_clean$core == "GeoB25202", , drop = FALSE]
geob_proxy$depth_in_core_cm <- as.numeric(geob_proxy$depth_in_core_cm)

geob_joined <- nearest_join_by_depth(
  reference_df = meta_geob_layer,
  reference_depth_col = "depth_in_core_cm",
  source_df = geob_proxy,
  source_depth_col = "depth_in_core_cm",
  tol_cm = 3
)
geob_joined <- geob_joined[geob_joined$proxy_match_within_tolerance, , drop = FALSE]

geob_proxy_vars <- c(
  "toc_percent_weight", "sum_expandable_clays", "sum_illites_micas", "sum_montmorillonites_smectites",
  "sum_greenland", "sum_iceland", "ratio_greenland_to_iceland", "d13c_per_mil"
)
geob_proxy_vars <- geob_proxy_vars[geob_proxy_vars %in% names(geob_joined)]

geob_rows <- vector("list", length(geob_proxy_vars))
for (i in seq_along(geob_proxy_vars)) {
  nm <- geob_proxy_vars[i]
  x <- as.numeric(geob_joined[[nm]])
  y <- as.numeric(geob_joined$library_concentration_log10)
  age <- as.numeric(geob_joined$y_bp)

  raw <- safe_spearman(x, y)

  ok <- is.finite(x) & is.finite(y) & is.finite(age)
  if (sum(ok) >= 12) {
    fit_y <- tryCatch(stats::lm(y[ok] ~ age[ok]), error = function(e) NULL)
    fit_x <- tryCatch(stats::lm(x[ok] ~ age[ok]), error = function(e) NULL)
    if (!is.null(fit_y) && !is.null(fit_x)) {
      detr <- safe_spearman(stats::residuals(fit_x), stats::residuals(fit_y), min_n = 12)
    } else {
      detr <- list(n = sum(ok), rho = NA_real_, p = NA_real_)
    }
  } else {
    detr <- list(n = sum(ok), rho = NA_real_, p = NA_real_)
  }

  geob_rows[[i]] <- data.frame(
    proxy = nm,
    n = raw$n,
    rho_raw = raw$rho,
    p_raw = raw$p,
    rho_age_adjusted_linear = detr$rho,
    p_age_adjusted_linear = detr$p,
    stringsAsFactors = FALSE
  )
}
geob_summary <- do.call(rbind, geob_rows)
geob_summary <- geob_summary[order(abs(geob_summary$rho_age_adjusted_linear), decreasing = TRUE), , drop = FALSE]
geob_summary$proxy_group <- ifelse(
  geob_summary$proxy %in% c("toc_percent_weight"), "TOC",
  ifelse(grepl("sum_|ratio_", geob_summary$proxy), "mineral/clay", "other")
)

if (nrow(geob_summary) > 0) {
  mineral_idx <- which(geob_summary$proxy_group == "mineral/clay")
  top_mineral <- if (length(mineral_idx) > 0) geob_summary$proxy[mineral_idx[1]] else NA_character_
  geob_summary$highlight_for_text <- geob_summary$proxy %in% c(top_mineral, "toc_percent_weight")
} else {
  geob_summary$highlight_for_text <- logical(0)
}

write_tsv(geob_summary, file.path(final_results, "geob25202_extra_proxy_summary.tsv"))

write_geob_proxy_plot(
  geob_summary,
  out_png = file.path(final_plots, "geob25202_extra_proxy_summary.png"),
  out_pdf = file.path(final_plots, "geob25202_extra_proxy_summary.pdf")
)

# ---- Curated figure set (selected only) ----
figure_map <- data.frame(
  figure_id = c(
    "supp_fig_st13_context",
    "supp_fig_st13_grouped_axes",
    "supp_fig_multicore_rbsr",
    "supp_fig_geob25202_extra"
  ),
  source = c(
    file.path("analysis", "st13_library_concentration_model", "plots", "st13_library_concentration_vs_age.png"),
    file.path("analysis", "st13_proxy_axes_model", "plots", "grouped_axis_coefficient_summary.png"),
    file.path("analysis", "multicore_rbsr_library_validation", "plots", "multicore_rbsr_comparison.png"),
    file.path(final_plots, "geob25202_extra_proxy_summary.png")
  ),
  dest = c(
    file.path(final_plots, "supp_fig_st13_context.png"),
    file.path(final_plots, "supp_fig_st13_grouped_axes.png"),
    file.path(final_plots, "supp_fig_multicore_rbsr.png"),
    file.path(final_plots, "supp_fig_geob25202_extra.png")
  ),
  title = c(
    "ST13 library concentration in age context",
    "ST13 grouped-axis coefficients",
    "Cross-core Rb/Sr comparison",
    "GeoB25202 extra-proxy summary"
  ),
  narrative_role = c(
    "Primary ST13 context",
    "Collinearity-resolved ST13 evidence",
    "External heterogeneity validation",
    "Convergent support within long core"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(figure_map))) {
  if (!file.exists(figure_map$source[i])) {
    warning("Figure source missing, skipped: ", figure_map$source[i])
    next
  }
  file.copy(figure_map$source[i], figure_map$dest[i], overwrite = TRUE)
}

figure_manifest <- figure_map[, c("figure_id", "dest", "title", "narrative_role")]
names(figure_manifest)[names(figure_manifest) == "dest"] <- "figure_file"
write_tsv(figure_manifest, file.path(final_results, "supplementary_figure_manifest.tsv"))

# ---- Curated table manifest ----
table_manifest <- data.frame(
  table_id = c(
    "supp_tab_exploratory_candidates",
    "supp_tab_st13_screening",
    "supp_tab_st13_primary_bootstrap",
    "supp_tab_st13_grouped_axes",
    "supp_tab_st13_state_dependence",
    "supp_tab_multicore_rbsr",
    "supp_tab_geob25202_extra"
  ),
  table_file = c(
    file.path(final_results, "exploratory_candidate_proxy_summary.tsv"),
    file.path(final_results, "st13_screening_selected_proxies.tsv"),
    file.path(final_results, "st13_primary_model_and_bootstrap_summary.tsv"),
    file.path(final_results, "st13_grouped_axis_model_summary.tsv"),
    file.path(final_results, "st13_state_dependence_interaction_summary.tsv"),
    file.path(final_results, "multicore_rbsr_validation_summary.tsv"),
    file.path(final_results, "geob25202_extra_proxy_summary.tsv")
  ),
  title = c(
    "Curated exploratory candidate proxies",
    "ST13 selected-proxy screening",
    "ST13 primary + block-bootstrap summary",
    "ST13 grouped-axis coefficients and fit",
    "ST13 state-dependence interaction summary",
    "Multicore Rb/Sr validation summary",
    "GeoB25202 extra-proxy effect summary"
  ),
  stringsAsFactors = FALSE
)
write_tsv(table_manifest, file.path(final_results, "supplementary_table_manifest.tsv"))

# ---- Key numbers for manuscript drafting ----
key_numbers <- data.frame(metric = character(0), value = character(0), stringsAsFactors = FALSE)
append_key <- function(metric, value) {
  key_numbers <<- rbind(key_numbers, data.frame(metric = metric, value = as.character(value), stringsAsFactors = FALSE))
}

top_rb_overall <- exploratory_candidates[exploratory_candidates$proxy_name == "rb_sr", , drop = FALSE]
if (nrow(top_rb_overall) > 0) {
  append_key("exploratory_rb_sr_rho", sprintf("%.3f", top_rb_overall$spearman_rho[1]))
  append_key("exploratory_rb_sr_q", signif(top_rb_overall$q_value_bh[1], 3))
}

st13_rb_primary <- st13_primary_boot[st13_primary_boot$model_id == "lm_primary" & st13_primary_boot$term == "z_rb_sr", , drop = FALSE]
if (nrow(st13_rb_primary) > 0) {
  append_key("st13_primary_z_rb_sr_estimate", sprintf("%.3f", st13_rb_primary$estimate[1]))
  append_key("st13_primary_z_rb_sr_p", signif(st13_rb_primary$p_value[1], 3))
}

st13_rb_boot <- st13_primary_boot[st13_primary_boot$model_id == "lm_block_bootstrap_primary" & st13_primary_boot$term == "z_rb_sr", , drop = FALSE]
if (nrow(st13_rb_boot) > 0) {
  append_key("st13_bootstrap_z_rb_sr_p", signif(st13_rb_boot$p_value[1], 3))
  append_key("st13_bootstrap_z_rb_sr_ci", sprintf("[%.3f, %.3f]", st13_rb_boot$ci_low[1], st13_rb_boot$ci_high[1]))
}

axis_terr <- axis_grouped[axis_grouped$term == "z_axis_terrigenous_mineralogical", , drop = FALSE]
if (nrow(axis_terr) > 0) {
  append_key("st13_grouped_axis_terrigenous_estimate", sprintf("%.3f", axis_terr$estimate[1]))
  append_key("st13_grouped_axis_terrigenous_p", signif(axis_terr$p_value[1], 3))
}

nest_terr <- axis_nested[axis_nested$model_id == "depth_plus_terrigenous", , drop = FALSE]
if (nrow(nest_terr) > 0) {
  append_key("st13_nested_terrigenous_model_p", signif(nest_terr$p_value[1], 3))
}

state_prod <- state_comparison[state_comparison$modifier_type == "productivity_carbonate", , drop = FALSE]
state_redox <- state_comparison[state_comparison$modifier_type == "redox", , drop = FALSE]
if (nrow(state_prod) > 0) append_key("st13_interaction_productivity_p", signif(state_prod$p_nested[1], 3))
if (nrow(state_redox) > 0) append_key("st13_interaction_redox_p", signif(state_redox$p_nested[1], 3))

for (cc in c("ST13", "GeoB25202")) {
  row <- multi_curated[multi_curated$core == cc, , drop = FALSE]
  if (nrow(row) > 0) {
    append_key(paste0("multicore_", cc, "_rb_sr_model_estimate"), sprintf("%.3f", row$rb_sr_model_estimate[1]))
    append_key(paste0("multicore_", cc, "_rb_sr_model_p"), signif(row$rb_sr_model_p[1], 3))
  }
}

if (nrow(geob_summary) > 0) {
  best <- geob_summary[1, , drop = FALSE]
  append_key("geob25202_strongest_extra_proxy", best$proxy[1])
  append_key("geob25202_strongest_extra_proxy_rho_age_adjusted", sprintf("%.3f", best$rho_age_adjusted_linear[1]))
  append_key("geob25202_strongest_extra_proxy_p_age_adjusted", signif(best$p_age_adjusted_linear[1], 3))
}

write_tsv(key_numbers, file.path(final_results, "supplementary_key_numbers.tsv"))

# ---- Supplementary narrative scaffold ----
outline_lines <- c(
  "Supplementary Results outline: library concentration vs proxy integration",
  "",
  "1) Exploratory candidate screening (minimal)",
  "   - Use: final/exploratory_candidate_proxy_summary.tsv",
  "   - Purpose: motivate ST13 terrigenous follow-up and GeoB25202 extra-proxy checks.",
  "",
  "2) ST13 primary evidence",
  "   - Use: final/st13_layer_table_essential.tsv",
  "   - Use: final/st13_screening_selected_proxies.tsv",
  "   - Use: final/st13_primary_model_and_bootstrap_summary.tsv",
  "   - Figure: plots/proxies/final/supp_fig_st13_context.png",
  "",
  "3) ST13 grouped-axis resolution of collinearity",
  "   - Use: final/st13_grouped_axis_membership_loadings.tsv",
  "   - Use: final/st13_grouped_axis_model_summary.tsv",
  "   - Use: final/st13_grouped_axis_nested_comparison.tsv",
  "   - Figure: plots/proxies/final/supp_fig_st13_grouped_axes.png",
  "",
  "4) ST13 state dependence (supporting negative/inconclusive)",
  "   - Use: final/st13_state_dependence_interaction_summary.tsv",
  "   - Interpretation: no robust interaction support from nested tests.",
  "",
  "5) Multicore external heterogeneity test (Rb/Sr)",
  "   - Use: final/multicore_rbsr_validation_summary.tsv",
  "   - Figure: plots/proxies/final/supp_fig_multicore_rbsr.png",
  "",
  "6) GeoB25202 convergent support from extra proxies",
  "   - Use: final/geob25202_extra_proxy_summary.tsv",
  "   - Figure: plots/proxies/final/supp_fig_geob25202_extra.png",
  "",
  "Cross-reference manuscript-ready numerical snippets in:",
  "   - final/supplementary_key_numbers.tsv"
)
writeLines(outline_lines, con = file.path(final_results, "supplementary_results_outline.txt"))

# ---- Assumptions, session, README ----
assumption_lines <- c(
  "Curated supplementary pipeline assumptions:",
  "1) Stage scripts in code/common are the authoritative source for model fitting logic.",
  "2) This script recomputes those stages and curates only manuscript-support outputs.",
  "3) Curated package is limited to results/common/proxies/final and plots/proxies/final.",
  "4) GeoB25202 extra-proxy effects are computed with nearest-depth matching (<=3 cm) and linear age detrending.",
  "5) Run from repository root so relative paths and provenance are deterministic."
)
writeLines(assumption_lines, con = file.path(final_results, "pipeline_assumptions.txt"))
writeLines(capture.output(sessionInfo()), con = file.path(final_results, "session_info.txt"))

readme_lines <- c(
  "# Curated supplementary package: library concentration vs proxies",
  "",
  "Generate with:",
  "",
  "```bash",
  "Rscript code/common/dna_concentration_exploration.R",
  "```",
  "",
  "## What this pipeline now does",
  "",
  "- Re-runs the necessary analysis stages (exploratory, ST13 primary, grouped axes, state dependence, multicore).",
  "- Curates a **small, supplement-facing** set of integrated tables/figures.",
  "- Avoids copying whole stage output trees into results/common/proxies/.",
  "",
  "## Final outputs used for Supplementary Results",
  "",
  "- `results/common/proxies/final/supplementary_table_manifest.tsv`",
  "- `results/common/proxies/final/supplementary_figure_manifest.tsv`",
  "- `results/common/proxies/final/supplementary_key_numbers.tsv`",
  "- `results/common/proxies/final/supplementary_results_outline.txt`",
  "- `plots/proxies/final/supp_fig_st13_context.png`",
  "- `plots/proxies/final/supp_fig_st13_grouped_axes.png`",
  "- `plots/proxies/final/supp_fig_multicore_rbsr.png`",
  "- `plots/proxies/final/supp_fig_geob25202_extra.png`",
  "",
  "## Provenance",
  "",
  "- Stage logs: `results/common/proxies/logs/`",
  "- Curated run manifest: `results/common/proxies/final/pipeline_run_manifest.tsv`",
  "- Assumptions: `results/common/proxies/final/pipeline_assumptions.txt`",
  "- Session info: `results/common/proxies/final/session_info.txt`",
  "",
  "## Guardrails",
  "",
  "- All results are observational and context/core specific.",
  "- State-dependence outputs are retained as concise negative/inconclusive support.",
  "- No downstream sequencing-output variables are used as explanatory proxies."
)
writeLines(readme_lines, con = file.path(results_root, "README.md"))

message("Curated supplementary pipeline complete.")
message("Final tables: ", final_results)
message("Final plots:  ", final_plots)
