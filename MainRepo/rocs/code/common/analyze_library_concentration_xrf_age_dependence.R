#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Paths and parameters
# -----------------------------------------------------------------------------
in_layer_table <- file.path(
  "results", "common", "proxies", "final",
  "library_concentration_xrf_supplementary_layer_table.tsv"
)

out_dir <- file.path("results", "common", "proxies", "final")
plot_dir <- file.path("plots", "proxies", "final")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_layer_table)) {
  stop("Missing required input: ", in_layer_table, "\nRun from repository root.")
}

params <- list(
  cores = c("ST13", "GeoB25202"),
  ratios = c("rb_sr", "ca_ti", "p_ti", "si_ti", "mn_fe", "al_si"),
  min_n_model = 10,
  min_n_cor = 8,
  interaction_alpha_hint = 0.10
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
safe_spearman <- function(x, y, min_n = 8) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < min_n || dplyr::n_distinct(x) < 3 || dplyr::n_distinct(y) < 3) {
    return(tibble(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(n = length(x), rho = unname(ct$estimate), p_value = ct$p.value)
}

safe_scale <- function(x) {
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  if (dplyr::n_distinct(x[is.finite(x)]) < 2) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

coef_row <- function(fit, term) {
  if (is.null(fit)) {
    return(tibble(
      estimate = NA_real_, std_error = NA_real_, p_value = NA_real_
    ))
  }
  tt <- summary(fit)$coefficients
  if (!(term %in% rownames(tt))) {
    return(tibble(
      estimate = NA_real_, std_error = NA_real_, p_value = NA_real_
    ))
  }
  tibble(
    estimate = tt[term, "Estimate"],
    std_error = tt[term, "Std. Error"],
    p_value = tt[term, "Pr(>|t|)"]
  )
}

fit_interaction_model <- function(df, ratio_name, core_name) {
  age_source <- if ("y_bp" %in% names(df) && sum(is.finite(df$y_bp)) >= params$min_n_model) "y_bp" else "depth_in_core_cm"

  dat <- df %>%
    select(core_clean, library_concentration_log10, depth_in_core_cm, all_of(age_source), all_of(ratio_name)) %>%
    filter(
      core_clean == core_name,
      is.finite(library_concentration_log10),
      is.finite(.data[[age_source]]),
      is.finite(.data[[ratio_name]])
    ) %>%
    mutate(
      age_value = .data[[age_source]],
      age_z = safe_scale(age_value),
      ratio_z = safe_scale(.data[[ratio_name]])
    ) %>%
    filter(is.finite(age_z), is.finite(ratio_z))

  if (nrow(dat) < params$min_n_model || dplyr::n_distinct(dat$ratio_z) < 4 || dplyr::n_distinct(dat$age_z) < 4) {
    return(tibble(
      core_clean = core_name,
      ratio_name = ratio_name,
      n_model = nrow(dat),
      age_source = age_source,
      age_term = "linear_z",
      interaction_term = "ratio_z:age_z",
      interaction_estimate = NA_real_,
      interaction_std_error = NA_real_,
      interaction_p_value = NA_real_,
      ratio_main_estimate = NA_real_,
      ratio_main_p_value = NA_real_,
      age_main_estimate = NA_real_,
      age_main_p_value = NA_real_,
      model_r_squared = NA_real_,
      note = "insufficient rows/variation for interaction model"
    ))
  }

  fit <- tryCatch(
    lm(library_concentration_log10 ~ ratio_z * age_z, data = dat),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      core_clean = core_name,
      ratio_name = ratio_name,
      n_model = nrow(dat),
      age_source = age_source,
      age_term = "linear_z",
      interaction_term = "ratio_z:age_z",
      interaction_estimate = NA_real_,
      interaction_std_error = NA_real_,
      interaction_p_value = NA_real_,
      ratio_main_estimate = NA_real_,
      ratio_main_p_value = NA_real_,
      age_main_estimate = NA_real_,
      age_main_p_value = NA_real_,
      model_r_squared = NA_real_,
      note = "lm fit failed"
    ))
  }

  int_term <- coef_row(fit, "ratio_z:age_z")
  ratio_term <- coef_row(fit, "ratio_z")
  age_coef <- coef_row(fit, "age_z")

  tibble(
    core_clean = core_name,
    ratio_name = ratio_name,
    n_model = nrow(dat),
    age_source = age_source,
    age_term = "linear_z",
    interaction_term = "ratio_z:age_z",
    interaction_estimate = int_term$estimate,
    interaction_std_error = int_term$std_error,
    interaction_p_value = int_term$p_value,
    ratio_main_estimate = ratio_term$estimate,
    ratio_main_p_value = ratio_term$p_value,
    age_main_estimate = age_coef$estimate,
    age_main_p_value = age_coef$p_value,
    model_r_squared = summary(fit)$r.squared,
    note = "interaction>0 implies stronger positive slope at older ages (or less negative)"
  )
}

make_strata_table <- function(df, ratio_name, core_name) {
  age_source <- if ("y_bp" %in% names(df) && sum(is.finite(df$y_bp)) >= params$min_n_model) "y_bp" else "depth_in_core_cm"

  dat <- df %>%
    filter(core_clean == core_name) %>%
    select(core_clean, library_concentration_log10, depth_in_core_cm, all_of(age_source), all_of(ratio_name)) %>%
    filter(
      is.finite(library_concentration_log10),
      is.finite(depth_in_core_cm),
      is.finite(.data[[age_source]]),
      is.finite(.data[[ratio_name]])
    ) %>%
    mutate(
      age_value = .data[[age_source]],
      half = ifelse(age_value <= median(age_value, na.rm = TRUE), "younger_half", "older_half"),
      tertile_id = ntile(age_value, 3),
      tertile = dplyr::case_when(
        tertile_id == 1 ~ "young_tertile",
        tertile_id == 2 ~ "mid_tertile",
        tertile_id == 3 ~ "old_tertile",
        TRUE ~ NA_character_
      )
    )

  if (nrow(dat) == 0) {
    return(tibble(
      core_clean = core_name,
      ratio_name = ratio_name,
      age_source = age_source,
      stratum_scheme = NA_character_,
      stratum = NA_character_,
      n = 0,
      spearman_raw = NA_real_,
      spearman_raw_p = NA_real_,
      depth_adjusted_rho = NA_real_,
      depth_adjusted_rho_p = NA_real_,
      depth_adjusted_slope = NA_real_,
      depth_adjusted_slope_p = NA_real_,
      note = "no rows"
    ))
  }

  strat_rows <- bind_rows(
    dat %>% transmute(stratum_scheme = "half", stratum = half, library_concentration_log10, depth_in_core_cm, ratio = .data[[ratio_name]]),
    dat %>% transmute(stratum_scheme = "tertile", stratum = tertile, library_concentration_log10, depth_in_core_cm, ratio = .data[[ratio_name]])
  )

  strat_rows %>%
    group_by(stratum_scheme, stratum) %>%
    group_modify(~ {
      d <- .x

      raw <- safe_spearman(d$ratio, d$library_concentration_log10, min_n = params$min_n_cor)

      depth_adj_rho <- tibble(n = nrow(d), rho = NA_real_, p_value = NA_real_)
      depth_adj_slope <- tibble(estimate = NA_real_, p_value = NA_real_)
      note <- "ok"

      if (
        nrow(d) >= params$min_n_cor &&
          dplyr::n_distinct(d$ratio) >= 3 &&
          dplyr::n_distinct(d$depth_in_core_cm) >= 3
      ) {
        y_res <- tryCatch(residuals(lm(library_concentration_log10 ~ scale(depth_in_core_cm), data = d)), error = function(e) rep(NA_real_, nrow(d)))
        x_res <- tryCatch(residuals(lm(ratio ~ scale(depth_in_core_cm), data = d)), error = function(e) rep(NA_real_, nrow(d)))
        depth_adj_rho <- safe_spearman(x_res, y_res, min_n = params$min_n_cor)

        d2 <- d %>% mutate(ratio_z = safe_scale(ratio), depth_z = safe_scale(depth_in_core_cm))
        fit <- tryCatch(
          lm(library_concentration_log10 ~ ratio_z + depth_z, data = d2),
          error = function(e) NULL
        )
        slope <- coef_row(fit, "ratio_z")
        depth_adj_slope <- tibble(estimate = slope$estimate, p_value = slope$p_value)
      } else {
        note <- "insufficient variation for depth-adjusted metrics"
      }

      tibble(
        n = nrow(d),
        spearman_raw = raw$rho,
        spearman_raw_p = raw$p_value,
        depth_adjusted_rho = depth_adj_rho$rho,
        depth_adjusted_rho_p = depth_adj_rho$p_value,
        depth_adjusted_slope = depth_adj_slope$estimate,
        depth_adjusted_slope_p = depth_adj_slope$p_value,
        note = note
      )
    }) %>%
    ungroup() %>%
    mutate(
      core_clean = core_name,
      ratio_name = ratio_name,
      age_source = age_source
    ) %>%
    select(core_clean, ratio_name, age_source, everything())
}

classify_signal <- function(interaction_p, old_slope, young_slope, old_rho, young_rho) {
  if (!is.finite(old_slope) || !is.finite(young_slope)) {
    return("no_clear_effect")
  }
  delta_abs <- abs(old_slope) - abs(young_slope)
  max_abs <- max(abs(old_slope), abs(young_slope), na.rm = TRUE)
  max_rho <- max(abs(old_rho), abs(young_rho), na.rm = TRUE)

  if (is.finite(interaction_p) && interaction_p <= params$interaction_alpha_hint) {
    if (delta_abs > 0.08) return("stronger_in_older_intervals")
    if (delta_abs < -0.08) return("stronger_in_younger_intervals")
    return("relatively_stable_effect")
  }

  if (max_abs < 0.12 && (!is.finite(max_rho) || max_rho < 0.20)) return("no_clear_effect")
  if (delta_abs > 0.15) return("stronger_in_older_intervals")
  if (delta_abs < -0.15) return("stronger_in_younger_intervals")
  if (abs(delta_abs) <= 0.08) return("relatively_stable_effect")
  "no_clear_effect"
}

make_core_plot <- function(df_long, core_name, outfile) {
  d <- df_long %>%
    filter(core_clean == core_name) %>%
    mutate(
      half = ifelse(age_value <= median(age_value, na.rm = TRUE), "younger_half", "older_half"),
      ratio_name = factor(ratio_name, levels = params$ratios)
    )

  p <- ggplot(d, aes(x = ratio_z, y = library_concentration_log10, color = half)) +
    geom_point(alpha = 0.65, size = 1.3) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    facet_wrap(~ratio_name, scales = "free_x", ncol = 2) +
    labs(
      title = paste0(core_name, ": age-stratified ratio signal vs log10(library concentration)"),
      subtitle = "Lines are younger-half vs older-half slopes; steeper older-half trend suggests age strengthening",
      x = "XRF ratio (z-score within core)",
      y = "log10(library concentration)",
      color = "Age stratum"
    ) +
    theme_bw()

  ggsave(outfile, p, width = 11, height = 8, dpi = 300)
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
layer <- read_tsv(in_layer_table, show_col_types = FALSE) %>%
  filter(core_clean %in% params$cores)

missing_ratio_cols <- params$ratios[!params$ratios %in% names(layer)]
if (length(missing_ratio_cols) > 0) {
  stop("Missing ratio columns in layer table: ", paste(missing_ratio_cols, collapse = ", "))
}

if (!"library_concentration_log10" %in% names(layer)) {
  stop("Expected column library_concentration_log10 in layer table.")
}

# -----------------------------------------------------------------------------
# A) Within-core age-interaction models
# -----------------------------------------------------------------------------
interaction_tbl <- purrr::map_dfr(params$cores, function(cc) {
  purrr::map_dfr(params$ratios, function(rr) fit_interaction_model(layer, rr, cc))
})

# -----------------------------------------------------------------------------
# B) Age-stratified analyses (half + tertiles)
# -----------------------------------------------------------------------------
strata_tbl <- purrr::map_dfr(params$cores, function(cc) {
  purrr::map_dfr(params$ratios, function(rr) make_strata_table(layer, rr, cc))
})

# -----------------------------------------------------------------------------
# C) Signal-strength summary table
# -----------------------------------------------------------------------------
half_tbl <- strata_tbl %>%
  filter(stratum_scheme == "half") %>%
  select(core_clean, ratio_name, stratum, n, spearman_raw, depth_adjusted_slope) %>%
  pivot_wider(
    names_from = stratum,
    values_from = c(n, spearman_raw, depth_adjusted_slope),
    names_sep = "__"
  )

tertile_tbl <- strata_tbl %>%
  filter(stratum_scheme == "tertile") %>%
  select(core_clean, ratio_name, stratum, depth_adjusted_slope) %>%
  pivot_wider(
    names_from = stratum,
    values_from = depth_adjusted_slope,
    names_prefix = "slope__"
  ) %>%
  mutate(
    tertile_abs_trend = dplyr::case_when(
      is.finite(slope__young_tertile) & is.finite(slope__mid_tertile) & is.finite(slope__old_tertile) &
        abs(slope__young_tertile) <= abs(slope__mid_tertile) & abs(slope__mid_tertile) <= abs(slope__old_tertile) ~ "increasing_abs_with_age",
      is.finite(slope__young_tertile) & is.finite(slope__mid_tertile) & is.finite(slope__old_tertile) &
        abs(slope__young_tertile) >= abs(slope__mid_tertile) & abs(slope__mid_tertile) >= abs(slope__old_tertile) ~ "decreasing_abs_with_age",
      is.finite(slope__young_tertile) & is.finite(slope__mid_tertile) & is.finite(slope__old_tertile) ~ "mixed_or_nonmonotonic",
      TRUE ~ "insufficient"
    )
  )

summary_tbl <- interaction_tbl %>%
  left_join(half_tbl, by = c("core_clean", "ratio_name")) %>%
  left_join(tertile_tbl, by = c("core_clean", "ratio_name")) %>%
  rowwise() %>%
  mutate(
    half_abs_delta_old_minus_young =
      abs(depth_adjusted_slope__older_half) - abs(depth_adjusted_slope__younger_half),
    signal_strength_call = classify_signal(
      interaction_p = interaction_p_value,
      old_slope = depth_adjusted_slope__older_half,
      young_slope = depth_adjusted_slope__younger_half,
      old_rho = spearman_raw__older_half,
      young_rho = spearman_raw__younger_half
    )
  ) %>%
  ungroup() %>%
  mutate(
    interaction_direction = dplyr::case_when(
      is.na(interaction_estimate) ~ NA_character_,
      interaction_estimate > 0 ~ "older_amplifies_positive_or_reduces_negative",
      interaction_estimate < 0 ~ "older_amplifies_negative_or_reduces_positive",
      TRUE ~ "neutral"
    )
  ) %>%
  select(
    core_clean, ratio_name,
    n_model, age_source,
    interaction_estimate, interaction_std_error, interaction_p_value,
    ratio_main_estimate, ratio_main_p_value,
    n__younger_half, n__older_half,
    spearman_raw__younger_half, spearman_raw__older_half,
    depth_adjusted_slope__younger_half, depth_adjusted_slope__older_half,
    half_abs_delta_old_minus_young,
    tertile_abs_trend,
    signal_strength_call,
    interaction_direction
  )

# -----------------------------------------------------------------------------
# D) Plots
# -----------------------------------------------------------------------------
plot_long <- layer %>%
  select(core_clean, library_concentration_log10, y_bp, all_of(params$ratios)) %>%
  pivot_longer(cols = all_of(params$ratios), names_to = "ratio_name", values_to = "ratio_value") %>%
  filter(
    core_clean %in% params$cores,
    is.finite(library_concentration_log10),
    is.finite(y_bp),
    is.finite(ratio_value)
  ) %>%
  group_by(core_clean, ratio_name) %>%
  mutate(
    ratio_z = safe_scale(ratio_value),
    age_value = y_bp
  ) %>%
  ungroup() %>%
  filter(is.finite(ratio_z))

make_core_plot(
  df_long = plot_long,
  core_name = "ST13",
  outfile = file.path(plot_dir, "library_concentration_xrf_age_dependence_st13.png")
)

make_core_plot(
  df_long = plot_long,
  core_name = "GeoB25202",
  outfile = file.path(plot_dir, "library_concentration_xrf_age_dependence_geob25202.png")
)

# -----------------------------------------------------------------------------
# E) Notes and outputs
# -----------------------------------------------------------------------------
notes_lines <- c(
  "Within-core age-effect modification analysis: XRF ratio signals vs library concentration",
  "",
  "Scope:",
  "- Cores analyzed independently: ST13 and GeoB25202.",
  "- Data source: supplementary layer table (results/common/proxies/final/library_concentration_xrf_supplementary_layer_table.tsv).",
  "- Ratios tested: rb_sr, ca_ti, p_ti, si_ti, mn_fe, al_si.",
  "",
  "Modeling approach:",
  "- For each core × ratio: lm(log10(library_concentration) ~ ratio_z * age_z).",
  "- Interaction coefficient (ratio_z:age_z) tests age dependence of the ratio slope.",
  "- Sensitivity stratification: younger vs older halves, plus young/mid/old tertiles.",
  "- Within each stratum, report raw Spearman rho and depth-adjusted standardized slope (lm with ratio_z + depth_z).",
  "",
  "Interpretation guardrails:",
  "- These are within-core observational effect-modification tests, not mechanistic proof.",
  "- Age and core span remain confounded with co-varying environmental and depositional changes.",
  "- A stronger old-section effect is consistent with (but does not prove) a shift from time-dominated to sediment-state-dominated control.",
  "",
  "Practical caveats:",
  "- GeoB25202 has sparse deep-age bins in places; interaction precision can be limited.",
  "- Depth-adjusted stratum metrics are simple local controls and do not remove all confounding."
)

write_tsv(
  interaction_tbl,
  file.path(out_dir, "library_concentration_xrf_age_dependence_interactions.tsv")
)

write_tsv(
  strata_tbl,
  file.path(out_dir, "library_concentration_xrf_age_dependence_strata.tsv")
)

write_tsv(
  summary_tbl,
  file.path(out_dir, "library_concentration_xrf_age_dependence_summary.tsv")
)

writeLines(
  notes_lines,
  con = file.path(out_dir, "library_concentration_xrf_age_dependence_notes.txt")
)

message("Age-dependence analysis complete.")
message("Outputs: ", out_dir)
message("Plots: ", plot_dir)
