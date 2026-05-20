#!/usr/bin/env Rscript

# Build a clean, paper-ready report package for Rb/Sr vs library_concentration.
# This script is intentionally provenance-first: it consolidates existing outputs
# from prior ST13 + multicore analyses rather than re-estimating new models.

options(stringsAsFactors = FALSE)

out_dir <- file.path("analysis", "rbsr_library_concentration_report")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

in_multicore <- file.path("analysis", "multicore_rbsr_library_validation")
in_st13 <- file.path("analysis", "st13_library_concentration_model")

required_files <- c(
  file.path(in_multicore, "multicore_rbsr_summary.tsv"),
  file.path(in_multicore, "multicore_rbsr_effect_details.tsv"),
  file.path(in_multicore, "per_core_age_stratified_rbsr_summary.tsv"),
  file.path(in_multicore, "per_core_rbsr_age_interaction_summary.tsv"),
  file.path(in_multicore, "per_core_temporal_vs_temporal_plus_rbsr_model_comparison.tsv"),
  file.path(in_multicore, "st13_layer_level_analysis_table.tsv"),
  file.path(in_multicore, "st5_layer_level_analysis_table.tsv"),
  file.path(in_multicore, "st8_layer_level_analysis_table.tsv"),
  file.path(in_multicore, "geob25202_layer_level_analysis_table.tsv"),
  file.path(in_st13, "valley_all_layers_with_residuals.tsv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input files. Run this script from repo root after the source analyses. Missing: ",
    paste(missing_files, collapse = ", ")
  )
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

safe_spearman <- function(x, y, min_n = 10) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < min_n || length(unique(x2)) < 3 || length(unique(y2)) < 3) {
    return(list(n = length(x2), rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(n = length(x2), rho = unname(ct$estimate), p = ct$p.value)
}

save_plot_dual <- function(stem, plot_fun, width_in = 9, height_in = 6, png_res = 300) {
  png_path <- file.path(plot_dir, paste0(stem, ".png"))
  pdf_path <- file.path(plot_dir, paste0(stem, ".pdf"))

  grDevices::png(
    filename = png_path,
    width = width_in,
    height = height_in,
    units = "in",
    res = png_res
  )
  plot_fun()
  grDevices::dev.off()

  grDevices::pdf(file = pdf_path, width = width_in, height = height_in, onefile = TRUE)
  plot_fun()
  grDevices::dev.off()

  c(png_path, pdf_path)
}

# -------------------------------
# Load source outputs (provenance inputs)
# -------------------------------
core_summary <- read_tsv(file.path(in_multicore, "multicore_rbsr_summary.tsv"))
effect_details <- read_tsv(file.path(in_multicore, "multicore_rbsr_effect_details.tsv"))
age_stratified <- read_tsv(file.path(in_multicore, "per_core_age_stratified_rbsr_summary.tsv"))
interaction_tbl <- read_tsv(file.path(in_multicore, "per_core_rbsr_age_interaction_summary.tsv"))
temporal_gain <- read_tsv(file.path(in_multicore, "per_core_temporal_vs_temporal_plus_rbsr_model_comparison.tsv"))

st13_layer <- read_tsv(file.path(in_multicore, "st13_layer_level_analysis_table.tsv"))
st5_layer <- read_tsv(file.path(in_multicore, "st5_layer_level_analysis_table.tsv"))
st8_layer <- read_tsv(file.path(in_multicore, "st8_layer_level_analysis_table.tsv"))
geob_layer <- read_tsv(file.path(in_multicore, "geob25202_layer_level_analysis_table.tsv"))
st13_valley <- read_tsv(file.path(in_st13, "valley_all_layers_with_residuals.tsv"))

# -------------------------------
# A) Summary tables
# -------------------------------
core_layers <- list(
  ST13 = st13_layer,
  ST5 = st5_layer,
  ST8 = st8_layer,
  GeoB25202 = geob_layer
)

core_age <- do.call(rbind, lapply(names(core_layers), function(core_nm) {
  d <- core_layers[[core_nm]]
  data.frame(
    core = core_nm,
    n_obs_age = sum(is.finite(d$y_bp)),
    age_min_y_bp = suppressWarnings(min(d$y_bp, na.rm = TRUE)),
    age_max_y_bp = suppressWarnings(max(d$y_bp, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}))

core_age$age_range_y_bp <- sprintf("%.1f to %.1f", core_age$age_min_y_bp, core_age$age_max_y_bp)

rbsr_effect_core <- effect_details[, c(
  "core",
  "rb_sr_detrended_p",
  "rb_sr_model_estimate",
  "rb_sr_model_p",
  "rb_sr_model_ci_low",
  "rb_sr_model_ci_high",
  "rb_sr_model_type"
), drop = FALSE]

rbsr_core_summary <- merge(core_summary, rbsr_effect_core, by = "core", all.x = TRUE, sort = FALSE)
rbsr_core_summary <- merge(rbsr_core_summary, core_age[, c("core", "age_min_y_bp", "age_max_y_bp", "age_range_y_bp")], by = "core", all.x = TRUE, sort = FALSE)

rbsr_core_summary <- rbsr_core_summary[, c(
  "core",
  "n_obs",
  "age_min_y_bp",
  "age_max_y_bp",
  "age_range_y_bp",
  "analysis_level",
  "rb_sr_raw_rho",
  "rb_sr_raw_p",
  "rb_sr_detrended_rho",
  "rb_sr_model_estimate",
  "rb_sr_model_p",
  "rb_sr_model_ci_low",
  "rb_sr_model_ci_high",
  "rb_sr_supported_yes_no",
  "notes"
)]

names(rbsr_core_summary) <- c(
  "core",
  "n_obs",
  "age_min_y_bp",
  "age_max_y_bp",
  "age_range_y_bp",
  "analysis_level",
  "rb_sr_raw_correlation",
  "rb_sr_raw_correlation_p",
  "rb_sr_age_adjusted_effect",
  "rb_sr_model_estimate",
  "rb_sr_model_p",
  "rb_sr_model_ci_low",
  "rb_sr_model_ci_high",
  "support_yes_no",
  "notes"
)

rbsr_core_summary$provenance_multicore_summary <- file.path(in_multicore, "multicore_rbsr_summary.tsv")
rbsr_core_summary$provenance_effect_details <- file.path(in_multicore, "multicore_rbsr_effect_details.tsv")

write_tsv(rbsr_core_summary, file.path(out_dir, "rbsr_core_summary.tsv"))

st13_strata <- age_stratified[
  age_stratified$core == "ST13" & age_stratified$scheme %in% c("median_half", "tertile"),
  , drop = FALSE
]

st13_tvp <- temporal_gain[
  temporal_gain$core == "ST13" & temporal_gain$scope == "age_stratum" &
    temporal_gain$scheme %in% c("median_half", "tertile"),
  c("scheme", "stratum", "delta_aic_temporal_minus_plus", "delta_corr_sq_temporal_minus_plus"),
  drop = FALSE
]

st13_age_strata_summary <- merge(st13_strata, st13_tvp, by = c("scheme", "stratum"), all.x = TRUE, sort = FALSE)

young_ref <- function(df, scheme_name) {
  s <- df[df$scheme == scheme_name, , drop = FALSE]
  if (scheme_name == "median_half") return(s$rb_sr_rho_age_adjusted_linear[s$stratum == "younger_half"][1])
  s$rb_sr_rho_age_adjusted_linear[s$stratum == "young_tertile"][1]
}

ref_med <- young_ref(st13_age_strata_summary, "median_half")
ref_ter <- young_ref(st13_age_strata_summary, "tertile")

st13_age_strata_summary$young_ref_age_adjusted_rho <- ifelse(
  st13_age_strata_summary$scheme == "median_half", ref_med, ref_ter
)
st13_age_strata_summary$delta_vs_young_ref_age_adjusted_rho <-
  st13_age_strata_summary$rb_sr_rho_age_adjusted_linear - st13_age_strata_summary$young_ref_age_adjusted_rho

st13_age_strata_summary$older_section_more_negative_than_young_ref <- NA_character_
is_old <- with(st13_age_strata_summary, (scheme == "median_half" & stratum == "older_half") |
  (scheme == "tertile" & stratum == "old_tertile"))
st13_age_strata_summary$older_section_more_negative_than_young_ref[is_old] <- ifelse(
  st13_age_strata_summary$delta_vs_young_ref_age_adjusted_rho[is_old] < 0,
  "yes", "no"
)

st13_age_strata_summary$support_note <- "stratum-specific estimate"
st13_age_strata_summary$support_note[st13_age_strata_summary$scheme == "tertile" & st13_age_strata_summary$stratum == "old_tertile"] <-
  "old tertile is more negative than young tertile (supportive older-section strengthening pattern)"
st13_age_strata_summary$support_note[st13_age_strata_summary$scheme == "median_half" & st13_age_strata_summary$stratum == "older_half"] <-
  "older half is not more negative than younger half in age-adjusted rho (pattern is mixed by stratification choice)"

st13_age_strata_summary <- st13_age_strata_summary[, c(
  "core", "trend_col", "scheme", "stratum", "n", "age_min", "age_max",
  "rb_sr_rho_raw", "rb_sr_p_raw",
  "rb_sr_rho_age_adjusted_linear", "rb_sr_p_age_adjusted_linear",
  "delta_aic_temporal_minus_plus", "delta_corr_sq_temporal_minus_plus",
  "young_ref_age_adjusted_rho", "delta_vs_young_ref_age_adjusted_rho",
  "older_section_more_negative_than_young_ref", "support_note"
)]

names(st13_age_strata_summary)[names(st13_age_strata_summary) == "rb_sr_rho_age_adjusted_linear"] <- "rb_sr_age_adjusted_effect"
names(st13_age_strata_summary)[names(st13_age_strata_summary) == "rb_sr_p_age_adjusted_linear"] <- "rb_sr_age_adjusted_effect_p"

st13_age_strata_summary$provenance_strata <- file.path(in_multicore, "per_core_age_stratified_rbsr_summary.tsv")
st13_age_strata_summary$provenance_model_gain <- file.path(in_multicore, "per_core_temporal_vs_temporal_plus_rbsr_model_comparison.tsv")

write_tsv(st13_age_strata_summary, file.path(out_dir, "rbsr_st13_age_strata_summary.tsv"))

whole_core_gain <- temporal_gain[
  temporal_gain$scope == "whole_core" & temporal_gain$stratum == "all",
  c("core", "delta_aic_temporal_minus_plus", "delta_corr_sq_temporal_minus_plus"),
  drop = FALSE
]

primary_core <- merge(
  rbsr_core_summary[, c(
    "core", "n_obs", "rb_sr_raw_correlation", "rb_sr_raw_correlation_p",
    "rb_sr_age_adjusted_effect", "rb_sr_model_estimate", "rb_sr_model_p",
    "rb_sr_model_ci_low", "rb_sr_model_ci_high", "support_yes_no"
  )],
  whole_core_gain,
  by = "core",
  all.x = TRUE,
  sort = FALSE
)
primary_core$scope <- "whole_core"
primary_core$scheme <- "none"
primary_core$stratum_or_comparison <- "all"
primary_core$interaction_beta <- NA_real_
primary_core$interaction_p <- NA_real_
primary_core$rb_sr_age_adjusted_effect_p <- NA_real_
primary_core$notes <- "Whole-core summary"

st13_interaction <- interaction_tbl[interaction_tbl$core == "ST13", , drop = FALSE]

st13_tertile <- st13_age_strata_summary[
  st13_age_strata_summary$scheme == "tertile",
  c(
    "n", "scheme", "stratum", "rb_sr_rho_raw", "rb_sr_p_raw", "rb_sr_age_adjusted_effect",
    "rb_sr_age_adjusted_effect_p", "delta_aic_temporal_minus_plus", "delta_corr_sq_temporal_minus_plus"
  ),
  drop = FALSE
]

if (nrow(st13_tertile) > 0) {
  st13_tertile$core <- "ST13"
  st13_tertile$n_obs <- st13_tertile$n
  st13_tertile$stratum_or_comparison <- st13_tertile$stratum
  st13_tertile$rb_sr_model_estimate <- NA_real_
  st13_tertile$rb_sr_model_p <- NA_real_
  st13_tertile$rb_sr_model_ci_low <- NA_real_
  st13_tertile$rb_sr_model_ci_high <- NA_real_
  st13_tertile$support_yes_no <- NA_character_
  st13_tertile$interaction_beta <- NA_real_
  st13_tertile$interaction_p <- NA_real_
  st13_tertile$notes <- "ST13 tertile stratified effects"
  st13_tertile <- st13_tertile[, c(
    "core", "n_obs", "rb_sr_rho_raw", "rb_sr_p_raw", "rb_sr_age_adjusted_effect",
    "rb_sr_age_adjusted_effect_p", "rb_sr_model_estimate", "rb_sr_model_p",
    "rb_sr_model_ci_low", "rb_sr_model_ci_high", "support_yes_no",
    "delta_aic_temporal_minus_plus", "delta_corr_sq_temporal_minus_plus",
    "scheme", "stratum_or_comparison", "interaction_beta", "interaction_p", "notes"
  )]
} else {
  st13_tertile <- data.frame()
}

if (nrow(st13_interaction) == 1) {
  inter_row <- data.frame(
    core = "ST13",
    n_obs = st13_interaction$n_obs,
    rb_sr_rho_raw = NA_real_,
    rb_sr_p_raw = NA_real_,
    rb_sr_age_adjusted_effect = NA_real_,
    rb_sr_age_adjusted_effect_p = NA_real_,
    rb_sr_model_estimate = st13_interaction$rbsr_main_beta,
    rb_sr_model_p = st13_interaction$rbsr_main_p,
    rb_sr_model_ci_low = NA_real_,
    rb_sr_model_ci_high = NA_real_,
    support_yes_no = NA_character_,
    delta_aic_temporal_minus_plus = NA_real_,
    delta_corr_sq_temporal_minus_plus = NA_real_,
    scheme = "interaction_model",
    stratum_or_comparison = "z_rb_sr:z_age_trend",
    interaction_beta = st13_interaction$rbsr_age_interaction_beta,
    interaction_p = st13_interaction$rbsr_age_interaction_p,
    notes = "ST13 interaction test; supportive if negative and small p, but here interpret as mixed evidence",
    stringsAsFactors = FALSE
  )
} else {
  inter_row <- data.frame()
}

primary_core$scheme <- "none"
primary_core$stratum_or_comparison <- "all"
primary_core$interaction_beta <- NA_real_
primary_core$interaction_p <- NA_real_

primary_core <- primary_core[, c(
  "core", "n_obs", "rb_sr_raw_correlation", "rb_sr_raw_correlation_p", "rb_sr_age_adjusted_effect",
  "rb_sr_age_adjusted_effect_p",
  "rb_sr_model_estimate", "rb_sr_model_p", "rb_sr_model_ci_low", "rb_sr_model_ci_high", "support_yes_no",
  "delta_aic_temporal_minus_plus", "delta_corr_sq_temporal_minus_plus",
  "scheme", "stratum_or_comparison", "interaction_beta", "interaction_p", "notes"
)]

names(primary_core)[names(primary_core) == "rb_sr_raw_correlation"] <- "rb_sr_rho_raw"
names(primary_core)[names(primary_core) == "rb_sr_raw_correlation_p"] <- "rb_sr_p_raw"

rbsr_primary_numbers_for_text <- rbind(
  primary_core,
  st13_tertile,
  inter_row
)

rbsr_primary_numbers_for_text$provenance <- ifelse(
  rbsr_primary_numbers_for_text$scheme %in% c("none", "interaction_model"),
  file.path(in_multicore, "multicore_rbsr_effect_details.tsv"),
  file.path(in_multicore, "per_core_age_stratified_rbsr_summary.tsv")
)

write_tsv(rbsr_primary_numbers_for_text, file.path(out_dir, "rbsr_primary_numbers_for_text.tsv"))

# -------------------------------
# B) Plots
# -------------------------------
save_plot_dual("cross_core_effect_summary", function() {
  d <- effect_details
  d <- d[order(d$rb_sr_model_estimate), , drop = FALSE]
  y <- seq_len(nrow(d))
  xlim <- range(c(d$rb_sr_model_ci_low, d$rb_sr_model_ci_high, d$rb_sr_detrended_rho, 0), na.rm = TRUE)
  if (!all(is.finite(xlim))) xlim <- c(-1, 1)

  plot(
    d$rb_sr_model_estimate, y,
    xlim = xlim,
    yaxt = "n",
    pch = 16,
    cex = 1.1,
    col = ifelse(d$rb_sr_model_estimate < 0, "#1f77b4", "#d62728"),
    xlab = "Age-adjusted/model-based Rb/Sr effect (standardized coefficient)",
    ylab = "",
    main = "Cross-core Rb/Sr effect on log10(library concentration)"
  )
  axis(2, at = y, labels = paste0(d$core, " (n=", d$n_obs, ")"), las = 1)
  segments(d$rb_sr_model_ci_low, y, d$rb_sr_model_ci_high, y, col = "#444444", lwd = 2)
  points(d$rb_sr_detrended_rho, y, pch = 17, col = "#2ca02c")
  abline(v = 0, lty = 2)
  legend(
    "bottomright",
    legend = c("Model estimate", "95% CI", "Detrended rho (reference)"),
    pch = c(16, NA, 17),
    lty = c(NA, 1, NA),
    col = c("#444444", "#444444", "#2ca02c"),
    bty = "n"
  )
})

save_plot_dual("st13_downcore_age_context", function() {
  st13 <- st13_layer[, c("layer_id", "y_bp", "library_concentration_log10"), drop = FALSE]
  st13 <- merge(st13, st13_valley[, c("layer_id", "valley_candidate")], by = "layer_id", all.x = TRUE, sort = FALSE)

  tert <- st13_age_strata_summary[st13_age_strata_summary$scheme == "tertile", , drop = FALSE]
  tert <- tert[order(tert$age_min), , drop = FALSE]

  xlim <- range(st13$y_bp, na.rm = TRUE)
  ylim <- range(st13$library_concentration_log10, na.rm = TRUE)

  plot(
    st13$y_bp, st13$library_concentration_log10,
    type = "n",
    xlab = "Age (y BP)",
    ylab = "log10(library concentration + pseudocount)",
    xlim = xlim,
    ylim = ylim,
    main = "ST13 age profile with valley context and age strata"
  )

  if (nrow(tert) > 0) {
    cols <- c(rgb(0.9, 0.9, 1, 0.4), rgb(0.9, 1, 0.9, 0.4), rgb(1, 0.9, 0.9, 0.4))
    for (i in seq_len(nrow(tert))) {
      rect(tert$age_min[i], ylim[1], tert$age_max[i], ylim[2], col = cols[i], border = NA)
    }
  }

  points(
    st13$y_bp, st13$library_concentration_log10,
    pch = ifelse(isTRUE(st13$valley_candidate), 17, 16),
    col = ifelse(isTRUE(st13$valley_candidate), "#d62728", "#1f77b4"),
    cex = 0.8
  )

  ord <- order(st13$y_bp)
  lw <- stats::lowess(st13$y_bp[ord], st13$library_concentration_log10[ord], f = 0.2)
  lines(lw, col = "black", lwd = 2)

  legend(
    "bottomright",
    legend = c("Non-valley layer", "Valley candidate", "Lowess trend", "Tertile background"),
    pch = c(16, 17, NA, 15),
    lty = c(NA, NA, 1, NA),
    col = c("#1f77b4", "#d62728", "black", "grey70"),
    bty = "n"
  )
})

save_plot_dual("st13_stratified_effect", function() {
  op <- par(mfrow = c(1, 2), mar = c(9, 4, 3, 1))

  s1 <- st13_age_strata_summary[st13_age_strata_summary$scheme == "median_half", , drop = FALSE]
  s1 <- s1[order(s1$stratum), , drop = FALSE]
  if (nrow(s1) > 0) {
    b1 <- barplot(
      s1$rb_sr_age_adjusted_effect,
      names.arg = s1$stratum,
      las = 2,
      col = c("#ff7f0e", "#1f77b4"),
      ylab = "Age-adjusted Spearman rho",
      main = "ST13 median-half strata"
    )
    abline(h = 0, lty = 2)
    text(b1, s1$rb_sr_age_adjusted_effect, labels = paste0("n=", s1$n), pos = ifelse(s1$rb_sr_age_adjusted_effect >= 0, 3, 1), cex = 0.8)
  } else {
    plot.new(); text(0.5, 0.5, "No median-half data")
  }

  s2 <- st13_age_strata_summary[st13_age_strata_summary$scheme == "tertile", , drop = FALSE]
  s2 <- s2[match(c("young_tertile", "mid_tertile", "old_tertile"), s2$stratum), , drop = FALSE]
  if (nrow(s2) > 0) {
    b2 <- barplot(
      s2$rb_sr_age_adjusted_effect,
      names.arg = s2$stratum,
      las = 2,
      col = c("#ffbb78", "#98df8a", "#1f77b4"),
      ylab = "Age-adjusted Spearman rho",
      main = "ST13 tertile strata"
    )
    abline(h = 0, lty = 2)
    text(b2, s2$rb_sr_age_adjusted_effect, labels = paste0("n=", s2$n), pos = ifelse(s2$rb_sr_age_adjusted_effect >= 0, 3, 1), cex = 0.8)
  } else {
    plot.new(); text(0.5, 0.5, "No tertile data")
  }

  par(op)
})

save_plot_dual("cross_core_raw_relationship_panel", function() {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  core_order <- c("ST13", "ST5", "ST8", "GeoB25202")
  for (core_nm in core_order) {
    d <- core_layers[[core_nm]]
    x <- d$rb_sr
    y <- d$library_concentration_log10
    ok <- is.finite(x) & is.finite(y)
    rr <- safe_spearman(x, y, min_n = 10)
    plot(
      x, y,
      pch = 16,
      col = "#1f77b4",
      cex = 0.7,
      xlab = "Rb/Sr",
      ylab = "log10(library concentration)",
      main = paste0(core_nm, " (n=", rr$n, ")")
    )
    if (sum(ok) >= 10) {
      abline(stats::lm(y[ok] ~ x[ok]), col = "#d62728", lwd = 1.8)
    }
    mtext(sprintf("Spearman rho = %.3f, p = %.2g", rr$rho, rr$p), side = 3, line = 0.1, cex = 0.8)
  }
  par(op)
})

save_plot_dual("cross_core_temporal_model_gain", function() {
  d <- temporal_gain[temporal_gain$scope == "whole_core" & temporal_gain$stratum == "all", , drop = FALSE]
  d <- d[match(c("ST13", "ST5", "ST8", "GeoB25202"), d$core), , drop = FALSE]
  op <- par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))

  if (nrow(d) > 0) {
    b1 <- barplot(
      d$delta_aic_temporal_minus_plus,
      names.arg = d$core,
      las = 2,
      col = "#1f77b4",
      ylab = "Delta AIC (temporal-only minus temporal+Rb/Sr)",
      main = "Whole-core model gain from adding Rb/Sr"
    )
    abline(h = 0, lty = 2)
    text(b1, d$delta_aic_temporal_minus_plus, labels = paste0("n=", d$n_obs), pos = ifelse(d$delta_aic_temporal_minus_plus >= 0, 3, 1), cex = 0.8)

    b2 <- barplot(
      d$delta_corr_sq_temporal_minus_plus,
      names.arg = d$core,
      las = 2,
      col = "#2ca02c",
      ylab = "Delta R2 (temporal+Rb/Sr minus temporal-only)",
      main = "Whole-core fit improvement"
    )
    abline(h = 0, lty = 2)
    text(b2, d$delta_corr_sq_temporal_minus_plus, labels = sprintf("%.3f", d$delta_corr_sq_temporal_minus_plus), pos = 3, cex = 0.8)
  } else {
    plot.new(); text(0.5, 0.5, "No whole-core temporal gain rows")
    plot.new(); text(0.5, 0.5, "No whole-core temporal gain rows")
  }

  par(op)
})

# -------------------------------
# C) Report text support files
# -------------------------------
get_num <- function(core_nm, col) {
  vv <- rbsr_core_summary[rbsr_core_summary$core == core_nm, col]
  if (length(vv) == 0) return(NA)
  vv[[1]]
}

st13_old_ter <- st13_age_strata_summary[
  st13_age_strata_summary$scheme == "tertile" & st13_age_strata_summary$stratum == "old_tertile",
  , drop = FALSE
]
st13_young_ter <- st13_age_strata_summary[
  st13_age_strata_summary$scheme == "tertile" & st13_age_strata_summary$stratum == "young_tertile",
  , drop = FALSE
]
st13_old_minus_young <- if (nrow(st13_old_ter) == 1 && nrow(st13_young_ter) == 1) {
  st13_old_ter$rb_sr_age_adjusted_effect[1] - st13_young_ter$rb_sr_age_adjusted_effect[1]
} else {
  NA_real_
}

supp_notes <- c(
  "Rb/Sr supplementary notes (consolidated from existing outputs):",
  "",
  sprintf("- ST13 whole-core: raw rho = %.3f (p=%.2g), age-adjusted rho = %.3f, model beta = %.3f (p=%.2g).",
          get_num("ST13", "rb_sr_raw_correlation"), get_num("ST13", "rb_sr_raw_correlation_p"),
          get_num("ST13", "rb_sr_age_adjusted_effect"), get_num("ST13", "rb_sr_model_estimate"), get_num("ST13", "rb_sr_model_p")),
  sprintf("- ST5 whole-core: raw rho = %.3f (p=%.2g), age-adjusted rho = %.3f, model beta = %.3f (p=%.2g).",
          get_num("ST5", "rb_sr_raw_correlation"), get_num("ST5", "rb_sr_raw_correlation_p"),
          get_num("ST5", "rb_sr_age_adjusted_effect"), get_num("ST5", "rb_sr_model_estimate"), get_num("ST5", "rb_sr_model_p")),
  sprintf("- ST8 whole-core: raw rho = %.3f (p=%.2g), age-adjusted rho = %.3f, model beta = %.3f (p=%.2g).",
          get_num("ST8", "rb_sr_raw_correlation"), get_num("ST8", "rb_sr_raw_correlation_p"),
          get_num("ST8", "rb_sr_age_adjusted_effect"), get_num("ST8", "rb_sr_model_estimate"), get_num("ST8", "rb_sr_model_p")),
  sprintf("- GeoB25202 whole-core: raw rho = %.3f (p=%.2g), age-adjusted rho = %.3f, model beta = %.3f (p=%.2g).",
          get_num("GeoB25202", "rb_sr_raw_correlation"), get_num("GeoB25202", "rb_sr_raw_correlation_p"),
          get_num("GeoB25202", "rb_sr_age_adjusted_effect"), get_num("GeoB25202", "rb_sr_model_estimate"), get_num("GeoB25202", "rb_sr_model_p")),
  "- Differential cross-core pattern: ST13 and GeoB25202 show negative supported effects; ST5 and ST8 do not show supported negative effects in the same framework.",
  sprintf("- ST13 tertiles: young age-adjusted rho = %.3f (p=%.2g), old age-adjusted rho = %.3f (p=%.2g), old-young difference = %.3f.",
          ifelse(nrow(st13_young_ter) == 1, st13_young_ter$rb_sr_age_adjusted_effect[1], NA_real_),
          ifelse(nrow(st13_young_ter) == 1, st13_young_ter$rb_sr_age_adjusted_effect_p[1], NA_real_),
          ifelse(nrow(st13_old_ter) == 1, st13_old_ter$rb_sr_age_adjusted_effect[1], NA_real_),
          ifelse(nrow(st13_old_ter) == 1, st13_old_ter$rb_sr_age_adjusted_effect_p[1], NA_real_),
          st13_old_minus_young),
  "- ST13 older-section strengthening is supportive in tertile stratification, but not uniformly supported across all stratification/model summaries; interpret as suggestive, not definitive.",
  "- Age has a strong structural role in all cores (temporal-only vs temporal+Rb/Sr comparisons included), but cross-core age span is confounded with core identity.",
  "- Guardrail: all results are observational and core-specific; these outputs support cautious interpretation, not mechanism proof.",
  "- Guardrail: GeoB25202 supports a strong negative Rb/Sr association but does not by itself demonstrate age-emergent strengthening."
)
writeLines(supp_notes, con = file.path(out_dir, "supplementary_results_notes.txt"))

methods_notes <- c(
  "Analysis methods notes (plain-language):",
  "",
  "- This package consolidates already-computed outputs from analysis/st13_library_concentration_model and analysis/multicore_rbsr_library_validation.",
  "- The target response is log10(library_concentration + pseudocount), as implemented in source analyses.",
  "- Core summaries combine: raw Rb/Sr correlation, age/depth-adjusted (detrended) Rb/Sr association, and depth/age-aware model coefficient when available.",
  "- ST13 age stratification summaries use existing younger/older halves and tertiles from per-core age-stratified outputs.",
  "- Cross-core comparison is descriptive/validation-oriented only; age range and core identity are confounded in multicore contrasts.",
  "- Downstream sequencing metrics are not used as explanatory variables in the source analyses and are not introduced here.",
  "- No new hypothesis tests or new model forms were introduced in this report script; this is a reproducible consolidation and figure packaging step."
)
writeLines(methods_notes, con = file.path(out_dir, "analysis_methods_notes.txt"))

readme_lines <- c(
  "# Rb/Sr library concentration report package",
  "",
  "This directory packages paper-ready tables/figures for the observed Rb/Sr effect on `library_concentration`, based on pre-existing analyses.",
  "",
  "## Provenance",
  "Primary upstream sources:",
  paste0("- `", file.path(in_multicore, "multicore_rbsr_summary.tsv"), "`"),
  paste0("- `", file.path(in_multicore, "multicore_rbsr_effect_details.tsv"), "`"),
  paste0("- `", file.path(in_multicore, "per_core_age_stratified_rbsr_summary.tsv"), "`"),
  paste0("- `", file.path(in_multicore, "per_core_rbsr_age_interaction_summary.tsv"), "`"),
  paste0("- `", file.path(in_multicore, "per_core_temporal_vs_temporal_plus_rbsr_model_comparison.tsv"), "`"),
  paste0("- `", file.path(in_st13, "valley_all_layers_with_residuals.tsv"), "`"),
  "",
  "## Tables",
  "- `rbsr_core_summary.tsv`: one-line-per-core summary with raw, age-adjusted, and model-based Rb/Sr effects + support flag.",
  "- `rbsr_st13_age_strata_summary.tsv`: ST13 younger/older half and tertile summaries, including older-vs-young contrast columns.",
  "- `rbsr_primary_numbers_for_text.tsv`: compact key numbers for supplementary prose drafting.",
  "",
  "## Plots (`plots/`, each in PNG + PDF)",
  "- `cross_core_effect_summary`: model-based and detrended Rb/Sr effect summary by core.",
  "- `st13_downcore_age_context`: ST13 library concentration vs age with valley candidates and age-strata context.",
  "- `st13_stratified_effect`: ST13 stratified age-adjusted Rb/Sr effects (median halves + tertiles).",
  "- `cross_core_raw_relationship_panel`: raw log-library vs Rb/Sr relationships by core.",
  "- `cross_core_temporal_model_gain`: temporal-only vs temporal+Rb/Sr gain across cores.",
  "",
  "## Text support files",
  "- `supplementary_results_notes.txt`: concise evidence-linked notes for manuscript supplement drafting.",
  "- `analysis_methods_notes.txt`: plain-language modeling/interpretation notes.",
  "",
  "## Rerun",
  "From repository root:",
  "",
  "```bash",
  "Rscript code/common/make_rbsr_library_concentration_report.R",
  "```",
  "",
  "## Guardrails carried through",
  "- Observational, core-specific associations only (not mechanism proof).",
  "- Cross-core age-span differences are confounded with core identity.",
  "- ST13 older-section strengthening is presented as supportive/mixed evidence, not definitive proof.",
  "- GeoB25202 supports a strong negative association, but not necessarily age-emergent strengthening."
)
writeLines(readme_lines, con = file.path(out_dir, "README.md"))

message("Rb/Sr report package complete: ", out_dir)
