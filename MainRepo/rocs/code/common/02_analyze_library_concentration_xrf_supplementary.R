#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(splines)
})

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Paths and analysis parameters
# -----------------------------------------------------------------------------
metadata_path <- file.path("data", "metadata_v5.tsv")
xrf_path <- file.path("data", "combined_xrf_not_normalized.tsv")

out_dir <- file.path("results", "common", "proxies", "final")
plot_dir <- file.path("plots", "proxies", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(metadata_path, xrf_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required inputs (run from repository root): ", paste(missing_files, collapse = ", "))
}

params <- list(
  min_n_cor = 8,
  min_n_model = 20,
  spline_df = 4,
  min_unique_trend = 6,
  core_targets = c("ST13", "ST5", "ST8", "GeoB25202"),
  interaction_ratios = c("rb_sr", "ca_ti")
)

ratio_definitions <- tibble::tribble(
  ~ratio_name, ~numerator, ~denominator, ~ratio_family, ~priority,
  "ba_ti", "ba", "ti", "productivity_carbonate", "required",
  "p_ti",  "p",  "ti", "productivity_carbonate", "required",
  "br_ti", "br", "ti", "productivity_carbonate", "required",
  "si_ti", "si", "ti", "productivity_carbonate", "required",
  "ca_ti", "ca", "ti", "productivity_carbonate", "required",
  "al_si", "al", "si", "terrigenous_mineralogical", "required",
  "ti_al", "ti", "al", "terrigenous_mineralogical", "required",
  "fe_al", "fe", "al", "terrigenous_mineralogical", "required",
  "rb_sr", "rb", "sr", "terrigenous_mineralogical", "required",
  "mn_fe", "mn", "fe", "redox", "required",
  "zr_al", "zr", "al", "terrigenous_mineralogical", "optional",
  "rb_al", "rb", "al", "terrigenous_mineralogical", "optional"
)

required_ratios <- ratio_definitions %>% filter(priority == "required") %>% pull(ratio_name)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
normalize_core <- function(x) stringr::str_replace(as.character(x), "_R[12]$", "")

safe_min_positive <- function(x, default = 1e-6) {
  vals <- x[is.finite(x) & x > 0]
  if (length(vals) == 0) return(default)
  min(vals)
}

choose_trend_var <- function(df) {
  if ("y_bp" %in% names(df)) {
    ok <- is.finite(df$y_bp)
    if (sum(ok) >= params$min_n_cor && dplyr::n_distinct(df$y_bp[ok]) >= params$min_unique_trend) return("y_bp")
  }
  if ("depth_in_core_cm" %in% names(df)) {
    ok <- is.finite(df$depth_in_core_cm)
    if (sum(ok) >= params$min_n_cor && dplyr::n_distinct(df$depth_in_core_cm[ok]) >= params$min_unique_trend) return("depth_in_core_cm")
  }
  NA_character_
}

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

safe_residualize <- function(y, trend, target_df = 4) {
  res <- rep(NA_real_, length(y))
  ok <- is.finite(y) & is.finite(trend)
  if (sum(ok) < params$min_n_cor) return(res)

  k <- dplyr::n_distinct(trend[ok])
  df_use <- max(1, min(target_df, k - 1))
  fit <- if (k >= 4) {
    tryCatch(lm(y[ok] ~ splines::ns(trend[ok], df = df_use)), error = function(e) NULL)
  } else {
    tryCatch(lm(y[ok] ~ trend[ok]), error = function(e) NULL)
  }
  if (is.null(fit)) return(res)
  res[ok] <- residuals(fit)
  res
}

derive_ratio <- function(df, num_col, den_col, ratio_col) {
  if (!(num_col %in% names(df)) || !(den_col %in% names(df))) {
    df[[ratio_col]] <- NA_real_
    return(df)
  }
  num <- suppressWarnings(as.numeric(df[[num_col]]))
  den <- suppressWarnings(as.numeric(df[[den_col]]))
  df[[ratio_col]] <- ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
  df
}

layer_intervals <- function(depths) {
  d <- sort(unique(depths[is.finite(depths)]))
  if (length(d) == 0) {
    return(tibble(depth_in_core_cm = numeric(), depth_lower = numeric(), depth_upper = numeric()))
  }
  if (length(d) == 1) {
    return(tibble(depth_in_core_cm = d, depth_lower = -Inf, depth_upper = Inf))
  }
  mids <- (d[-1] + d[-length(d)]) / 2
  tibble(
    depth_in_core_cm = d,
    depth_lower = c(-Inf, mids),
    depth_upper = c(mids, Inf)
  )
}

aggregate_xrf_to_layers <- function(layer_df, xrf_core_df, ratio_names) {
  ints <- layer_intervals(layer_df$depth_in_core_cm)
  if (nrow(ints) == 0) return(tibble())

  out <- ints %>% select(depth_in_core_cm, depth_lower, depth_upper)
  out$n_xrf_rows_interval <- 0L
  for (rn in ratio_names) {
    out[[rn]] <- NA_real_
    out[[paste0("n_valid_", rn)]] <- NA_integer_
  }

  for (i in seq_len(nrow(out))) {
    lo <- out$depth_lower[i]
    hi <- out$depth_upper[i]
    idx <- which(is.finite(xrf_core_df$depth_in_core_cm) & xrf_core_df$depth_in_core_cm > lo & xrf_core_df$depth_in_core_cm <= hi)
    out$n_xrf_rows_interval[i] <- length(idx)
    if (length(idx) == 0) next

    for (rn in ratio_names) {
      vals <- xrf_core_df[[rn]][idx]
      vals_ok <- vals[is.finite(vals)]
      out[[paste0("n_valid_", rn)]][i] <- length(vals_ok)
      if (length(vals_ok) > 0) out[[rn]][i] <- median(vals_ok)
    }
  }

  out
}

fit_ratio_model <- function(df, ratio_name, core_adjusted = FALSE, try_gls = FALSE) {
  trend_var <- choose_trend_var(df)
  if (is.na(trend_var)) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = NA_character_,
      model_type = "not_fit",
      n_model = 0,
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "no usable trend variable"
    ))
  }

  cols <- c("library_concentration_log10", trend_var, ratio_name)
  if (core_adjusted) cols <- c(cols, "core_clean")

  dat <- df %>%
    select(all_of(cols)) %>%
    filter(is.finite(library_concentration_log10), is.finite(.data[[trend_var]]), is.finite(.data[[ratio_name]])) %>%
    mutate(ratio_z = as.numeric(scale(.data[[ratio_name]])))

  if (nrow(dat) < params$min_n_model || dplyr::n_distinct(dat$ratio_z) < 4 || !all(is.finite(dat$ratio_z))) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_var,
      model_type = "not_fit",
      n_model = nrow(dat),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "insufficient rows or ratio variation"
    ))
  }

  k <- dplyr::n_distinct(dat[[trend_var]])
  df_ns <- max(1, min(params$spline_df, k - 1))
  rhs <- paste0("splines::ns(", trend_var, ", df = ", df_ns, ") + ratio_z")
  if (core_adjusted) rhs <- paste(rhs, "+ core_clean")
  fml <- as.formula(paste("library_concentration_log10 ~", rhs))

  if (try_gls && requireNamespace("nlme", quietly = TRUE)) {
    dat_gls <- dat %>% arrange(.data[[trend_var]])
    gls_fit <- tryCatch(
      nlme::gls(
        model = fml,
        data = dat_gls,
        method = "REML",
        correlation = nlme::corCAR1(form = stats::as.formula(paste0("~", trend_var)))
      ),
      error = function(e) NULL
    )
    if (!is.null(gls_fit)) {
      tt <- summary(gls_fit)$tTable
      if ("ratio_z" %in% rownames(tt)) {
        ci <- tryCatch(nlme::intervals(gls_fit, which = "coef")$coef, error = function(e) NULL)
        return(tibble(
          ratio_name = ratio_name,
          trend_used = trend_var,
          model_type = "gls_corCAR1",
          n_model = nrow(dat_gls),
          estimate = tt["ratio_z", "Value"],
          std_error = tt["ratio_z", "Std.Error"],
          p_value = tt["ratio_z", "p-value"],
          ci_low = if (!is.null(ci)) ci["ratio_z", "lower"] else NA_real_,
          ci_high = if (!is.null(ci)) ci["ratio_z", "upper"] else NA_real_,
          note = paste0("trend spline df=", df_ns)
        ))
      }
    }
  }

  lm_fit <- tryCatch(lm(fml, data = dat), error = function(e) NULL)
  if (is.null(lm_fit)) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_var,
      model_type = "not_fit",
      n_model = nrow(dat),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "lm fit failed"
    ))
  }

  tt <- summary(lm_fit)$coefficients
  if (!("ratio_z" %in% rownames(tt))) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_var,
      model_type = "lm",
      n_model = nrow(dat),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "ratio term not estimable"
    ))
  }

  ci <- tryCatch(confint(lm_fit), error = function(e) NULL)
  tibble(
    ratio_name = ratio_name,
    trend_used = trend_var,
    model_type = "lm",
    n_model = nrow(dat),
    estimate = tt["ratio_z", "Estimate"],
    std_error = tt["ratio_z", "Std. Error"],
    p_value = tt["ratio_z", "Pr(>|t|)"],
    ci_low = if (!is.null(ci)) ci["ratio_z", 1] else NA_real_,
    ci_high = if (!is.null(ci)) ci["ratio_z", 2] else NA_real_,
    note = paste0("trend spline df=", df_ns)
  )
}

compute_interaction_slopes <- function(df, ratio_name) {
  trend_var <- choose_trend_var(df)
  if (is.na(trend_var)) {
    return(list(summary = tibble(), detail = tibble()))
  }

  dat <- df %>%
    select(core_clean, library_concentration_log10, all_of(trend_var), all_of(ratio_name)) %>%
    filter(is.finite(library_concentration_log10), is.finite(.data[[trend_var]]), is.finite(.data[[ratio_name]])) %>%
    mutate(
      core_clean = factor(core_clean),
      ratio_z = as.numeric(scale(.data[[ratio_name]]))
    )

  if (nrow(dat) < params$min_n_model || dplyr::n_distinct(dat$core_clean) < 2 || !all(is.finite(dat$ratio_z))) {
    return(list(summary = tibble(), detail = tibble()))
  }

  k <- dplyr::n_distinct(dat[[trend_var]])
  df_ns <- max(1, min(params$spline_df, k - 1))

  fml0 <- as.formula(paste0("library_concentration_log10 ~ core_clean + splines::ns(", trend_var, ", df = ", df_ns, ") + ratio_z"))
  fml1 <- as.formula(paste0("library_concentration_log10 ~ core_clean + splines::ns(", trend_var, ", df = ", df_ns, ") + ratio_z * core_clean"))

  m0 <- tryCatch(lm(fml0, data = dat), error = function(e) NULL)
  m1 <- tryCatch(lm(fml1, data = dat), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(list(summary = tibble(), detail = tibble()))

  cmp <- tryCatch(anova(m0, m1), error = function(e) NULL)
  interaction_p <- if (!is.null(cmp) && nrow(cmp) >= 2) cmp$`Pr(>F)`[2] else NA_real_

  coef_vec <- coef(m1)
  V <- vcov(m1)
  cn <- names(coef_vec)
  baseline_core <- levels(dat$core_clean)[1]

  core_slopes <- map_dfr(levels(dat$core_clean), function(cc) {
    g <- rep(0, length(coef_vec))
    names(g) <- cn
    if ("ratio_z" %in% cn) g["ratio_z"] <- 1

    inter_name <- paste0("ratio_z:core_clean", cc)
    inter_name_rev <- paste0("core_clean", cc, ":ratio_z")
    if (inter_name %in% cn) g[inter_name] <- 1
    if (inter_name_rev %in% cn) g[inter_name_rev] <- 1

    est <- as.numeric(sum(g * coef_vec))
    se <- sqrt(as.numeric(t(g) %*% V %*% g))
    z <- ifelse(is.finite(se) && se > 0, est / se, NA_real_)
    p <- ifelse(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_)
    ci_low <- est - 1.96 * se
    ci_high <- est + 1.96 * se

    tibble(
      ratio_name = ratio_name,
      core_clean = cc,
      trend_used = trend_var,
      n_model = nrow(dat),
      interaction_global_p = interaction_p,
      estimate = est,
      std_error = se,
      p_value = p,
      ci_low = ci_low,
      ci_high = ci_high,
      baseline_core = baseline_core
    )
  })

  summary_tbl <- core_slopes %>%
    summarise(
      ratio_name = first(ratio_name),
      trend_used = first(trend_used),
      n_model = first(n_model),
      interaction_global_p = first(interaction_global_p),
      n_cores = n(),
      n_positive_slopes = sum(is.finite(estimate) & estimate > 0),
      n_negative_slopes = sum(is.finite(estimate) & estimate < 0),
      heterogeneity_call = dplyr::case_when(
        is.na(interaction_global_p) ~ "inconclusive",
        interaction_global_p <= 0.05 ~ "strong_heterogeneity",
        interaction_global_p <= 0.10 ~ "moderate_heterogeneity",
        TRUE ~ "limited_heterogeneity"
      )
    )

  list(summary = summary_tbl, detail = core_slopes)
}

save_plot_both <- function(plot_obj, filename_base, width = 10, height = 6) {
  ggsave(file.path(plot_dir, paste0(filename_base, ".png")), plot_obj, width = width, height = height, dpi = 300)
  ggsave(file.path(plot_dir, paste0(filename_base, ".pdf")), plot_obj, width = width, height = height)
}

# -----------------------------------------------------------------------------
# 1) Load and collapse metadata to one row per sediment layer
# -----------------------------------------------------------------------------
metadata_raw <- read_tsv(metadata_path, show_col_types = FALSE)
xrf_raw <- read_tsv(xrf_path, show_col_types = FALSE)

metadata_layers <- metadata_raw %>%
  mutate(core_clean = normalize_core(core)) %>%
  filter(core_clean %in% params$core_targets, is.finite(depth_in_core_cm), is.finite(y_bp)) %>%
  group_by(core_clean, depth_in_core_cm, y_bp) %>%
  summarise(
    n_technical_replicates = n(),
    n_library_nonmissing = sum(is.finite(library_concentration)),
    library_concentration_mean = mean(library_concentration, na.rm = TRUE),
    library_concentration_median = median(library_concentration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    library_concentration_mean = ifelse(is.nan(library_concentration_mean), NA_real_, library_concentration_mean),
    library_concentration_median = ifelse(is.nan(library_concentration_median), NA_real_, library_concentration_median)
  )

pseudocount <- safe_min_positive(metadata_layers$library_concentration_mean, default = 1e-6) / 2

metadata_layers <- metadata_layers %>%
  mutate(
    library_concentration_log10 = log10(library_concentration_mean + pseudocount),
    layer_id = paste(core_clean, depth_in_core_cm, y_bp, sep = "|")
  ) %>%
  arrange(core_clean, depth_in_core_cm)

# -----------------------------------------------------------------------------
# 2) Derive ratio-only XRF features from non-normalized source table
# -----------------------------------------------------------------------------
xrf <- xrf_raw %>%
  mutate(core_clean = as.character(core))

for (i in seq_len(nrow(ratio_definitions))) {
  xrf <- derive_ratio(
    xrf,
    num_col = ratio_definitions$numerator[i],
    den_col = ratio_definitions$denominator[i],
    ratio_col = ratio_definitions$ratio_name[i]
  )
}

ratio_names <- ratio_definitions$ratio_name

# -----------------------------------------------------------------------------
# 3) Aggregate dense XRF to library-layer intervals (median)
# -----------------------------------------------------------------------------
layer_table <- map_dfr(params$core_targets, function(cc) {
  layers <- metadata_layers %>% filter(core_clean == cc)
  xrf_core <- xrf %>% filter(core_clean == cc, is.finite(depth_in_core_cm))
  if (nrow(layers) == 0) return(tibble())

  agg <- aggregate_xrf_to_layers(layers, xrf_core, ratio_names)
  layers %>%
    left_join(agg, by = "depth_in_core_cm")
}) %>%
  mutate(core_clean = as.character(core_clean))

# -----------------------------------------------------------------------------
# 4) QC summary and coverage diagnostics
# -----------------------------------------------------------------------------
qc_summary <- map_dfr(params$core_targets, function(cc) {
  core_layers <- layer_table %>% filter(core_clean == cc)
  core_xrf <- xrf %>% filter(core_clean == cc)

  map_dfr(seq_len(nrow(ratio_definitions)), function(i) {
    rn <- ratio_definitions$ratio_name[i]
    num <- ratio_definitions$numerator[i]
    den <- ratio_definitions$denominator[i]
    num_vals <- if (num %in% names(core_xrf)) suppressWarnings(as.numeric(core_xrf[[num]])) else rep(NA_real_, nrow(core_xrf))
    den_vals <- if (den %in% names(core_xrf)) suppressWarnings(as.numeric(core_xrf[[den]])) else rep(NA_real_, nrow(core_xrf))

    tibble(
      core_clean = cc,
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[i],
      priority = ratio_definitions$priority[i],
      n_layers_core = nrow(core_layers),
      n_xrf_rows_core = nrow(core_xrf),
      n_layers_with_response = sum(is.finite(core_layers$library_concentration_log10)),
      n_layers_with_ratio = sum(is.finite(core_layers[[rn]])),
      pct_layers_with_ratio = 100 * mean(is.finite(core_layers[[rn]])),
      median_xrf_rows_per_interval = median(core_layers$n_xrf_rows_interval, na.rm = TRUE),
      median_valid_ratio_per_interval = median(core_layers[[paste0("n_valid_", rn)]], na.rm = TRUE),
      denominator_nonpositive_or_missing = sum(!is.finite(den_vals) | den_vals <= 0),
      numerator_missing = sum(!is.finite(num_vals))
    )
  })
}) %>%
  mutate(
    median_xrf_rows_per_interval = ifelse(is.nan(median_xrf_rows_per_interval), NA_real_, median_xrf_rows_per_interval),
    median_valid_ratio_per_interval = ifelse(is.nan(median_valid_ratio_per_interval), NA_real_, median_valid_ratio_per_interval)
  ) %>%
  arrange(core_clean, ratio_family, ratio_name)

# -----------------------------------------------------------------------------
# 5) Per-core associations (raw and detrended Spearman)
# -----------------------------------------------------------------------------
core_associations <- map_dfr(params$core_targets, function(cc) {
  dat <- layer_table %>% filter(core_clean == cc)
  trend_var <- choose_trend_var(dat)

  map_dfr(seq_len(nrow(ratio_definitions)), function(i) {
    rn <- ratio_definitions$ratio_name[i]
    raw <- safe_spearman(dat[[rn]], dat$library_concentration_log10, min_n = params$min_n_cor)

    det <- tibble(n = NA_integer_, rho = NA_real_, p_value = NA_real_)
    if (!is.na(trend_var)) {
      xr <- safe_residualize(dat[[rn]], dat[[trend_var]], target_df = params$spline_df)
      yr <- safe_residualize(dat$library_concentration_log10, dat[[trend_var]], target_df = params$spline_df)
      det <- safe_spearman(xr, yr, min_n = params$min_n_cor)
    }

    tibble(
      core_clean = cc,
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[i],
      priority = ratio_definitions$priority[i],
      trend_used = trend_var,
      n_raw = raw$n,
      spearman_raw = raw$rho,
      p_raw = raw$p_value,
      n_detrended = det$n,
      spearman_detrended = det$rho,
      p_detrended = det$p_value
    )
  })
}) %>%
  group_by(core_clean) %>%
  mutate(
    q_raw_bh_core = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh_core = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 6) Per-core depth/age adjusted models (required ratios)
# -----------------------------------------------------------------------------
core_models <- map_dfr(params$core_targets, function(cc) {
  dat <- layer_table %>% filter(core_clean == cc)
  gls_ok <- cc %in% c("ST13", "GeoB25202")

  map_dfr(required_ratios, function(rn) {
    fit_ratio_model(dat, ratio_name = rn, core_adjusted = FALSE, try_gls = gls_ok) %>%
      mutate(
        core_clean = cc,
        ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
        priority = ratio_definitions$priority[match(rn, ratio_definitions$ratio_name)],
        estimate_direction = case_when(
          is.na(estimate) ~ NA_character_,
          estimate > 0 ~ "positive",
          estimate < 0 ~ "negative",
          TRUE ~ "zero"
        )
      )
  })
}) %>%
  group_by(core_clean) %>%
  mutate(q_model_bh_core = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH"))) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 7) Cross-core pooled analysis
# -----------------------------------------------------------------------------
pooled_trend_var <- choose_trend_var(layer_table)

pooled_models <- map_dfr(required_ratios, function(rn) {
  fit <- fit_ratio_model(layer_table, ratio_name = rn, core_adjusted = TRUE, try_gls = FALSE)
  fit %>%
    mutate(
      ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
      priority = ratio_definitions$priority[match(rn, ratio_definitions$ratio_name)],
      analysis_type = "pooled_core_plus_smooth_trend_lm"
    )
}) %>%
  mutate(
    estimate_direction = case_when(
      is.na(estimate) ~ NA_character_,
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "zero"
    )
  )

pooled_residual_correlations <- map_dfr(required_ratios, function(rn) {
  dat <- layer_table %>%
    select(core_clean, library_concentration_log10, y_bp, depth_in_core_cm, all_of(rn)) %>%
    filter(is.finite(library_concentration_log10), is.finite(.data[[rn]]))

  trend_var <- choose_trend_var(dat)
  if (nrow(dat) < params$min_n_cor || is.na(trend_var)) {
    return(tibble(
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
      trend_used = trend_var,
      n_raw = nrow(dat),
      spearman_raw = NA_real_,
      p_raw = NA_real_,
      n_resid = nrow(dat),
      spearman_resid = NA_real_,
      p_resid = NA_real_
    ))
  }

  k <- dplyr::n_distinct(dat[[trend_var]])
  df_ns <- max(1, min(params$spline_df, k - 1))

  fy <- as.formula(paste0("library_concentration_log10 ~ core_clean + splines::ns(", trend_var, ", df = ", df_ns, ")"))
  fx <- as.formula(paste0(rn, " ~ core_clean + splines::ns(", trend_var, ", df = ", df_ns, ")"))

  my <- tryCatch(lm(fy, data = dat), error = function(e) NULL)
  mx <- tryCatch(lm(fx, data = dat), error = function(e) NULL)

  raw <- safe_spearman(dat[[rn]], dat$library_concentration_log10, min_n = params$min_n_cor)
  if (is.null(my) || is.null(mx)) {
    return(tibble(
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
      trend_used = trend_var,
      n_raw = raw$n,
      spearman_raw = raw$rho,
      p_raw = raw$p_value,
      n_resid = nrow(dat),
      spearman_resid = NA_real_,
      p_resid = NA_real_
    ))
  }

  resid_cor <- safe_spearman(residuals(mx), residuals(my), min_n = params$min_n_cor)
  tibble(
    ratio_name = rn,
    ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
    trend_used = trend_var,
    n_raw = raw$n,
    spearman_raw = raw$rho,
    p_raw = raw$p_value,
    n_resid = resid_cor$n,
    spearman_resid = resid_cor$rho,
    p_resid = resid_cor$p_value
  )
}) %>%
  mutate(
    q_resid_bh = ifelse(is.na(p_resid), NA_real_, p.adjust(p_resid, method = "BH"))
  )

interaction_results <- map(params$interaction_ratios, ~ compute_interaction_slopes(layer_table, .x))
interaction_summary <- bind_rows(map(interaction_results, "summary"))
interaction_detail <- bind_rows(map(interaction_results, "detail")) %>%
  left_join(ratio_definitions %>% select(ratio_name, ratio_family), by = "ratio_name")

if (nrow(interaction_detail) > 0) {
  interaction_detail <- interaction_detail %>%
    group_by(ratio_name) %>%
    mutate(q_core_slope_bh = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH"))) %>%
    ungroup()
}

pooled_models <- pooled_models %>%
  left_join(
    interaction_summary %>% select(ratio_name, interaction_global_p, heterogeneity_call),
    by = "ratio_name"
  ) %>%
  mutate(
    consistency_call = case_when(
      is.na(estimate) ~ "inconclusive",
      !is.na(interaction_global_p) & interaction_global_p <= 0.10 ~ "core_specific_or_heterogeneous",
      estimate > 0 ~ "consistent_positive_direction",
      estimate < 0 ~ "consistent_negative_direction",
      TRUE ~ "weak_or_null"
    )
  )

# -----------------------------------------------------------------------------
# 8) Ratio-family synthesis
# -----------------------------------------------------------------------------
family_core_direction <- core_models %>%
  mutate(sign_num = case_when(estimate > 0 ~ 1, estimate < 0 ~ -1, TRUE ~ 0)) %>%
  group_by(ratio_family) %>%
  summarise(
    n_core_coefficients = sum(is.finite(estimate)),
    n_positive = sum(sign_num > 0, na.rm = TRUE),
    n_negative = sum(sign_num < 0, na.rm = TRUE),
    n_nominal_support = sum(is.finite(p_value) & p_value <= 0.10),
    .groups = "drop"
  )

family_pooled <- pooled_models %>%
  group_by(ratio_family) %>%
  summarise(
    n_ratios = n(),
    mean_pooled_estimate = mean(estimate, na.rm = TRUE),
    median_pooled_estimate = median(estimate, na.rm = TRUE),
    n_positive_pooled = sum(is.finite(estimate) & estimate > 0),
    n_negative_pooled = sum(is.finite(estimate) & estimate < 0),
    .groups = "drop"
  )

ratio_family_summary <- family_pooled %>%
  left_join(family_core_direction, by = "ratio_family") %>%
  mutate(
    family_direction_call = case_when(
      is.nan(mean_pooled_estimate) ~ "inconclusive",
      mean_pooled_estimate > 0 ~ "net_positive",
      mean_pooled_estimate < 0 ~ "net_negative",
      TRUE ~ "neutral"
    ),
    sign_consistency_index = ifelse(
      (n_positive + n_negative) > 0,
      pmax(n_positive, n_negative) / (n_positive + n_negative),
      NA_real_
    )
  )

# -----------------------------------------------------------------------------
# 9) Key numbers and narrative notes
# -----------------------------------------------------------------------------
n_layers_total <- nrow(layer_table)
n_layers_by_core <- layer_table %>% count(core_clean, name = "n_layers")

key_ratios <- c("ca_ti", "rb_sr", "mn_fe", "al_si", "ti_al", "fe_al")
key_pooled <- pooled_models %>% filter(ratio_name %in% key_ratios)

key_numbers <- bind_rows(
  tibble(metric = "n_layers_total", value = n_layers_total, note = "Layer-level rows after replicate collapse and XRF aggregation"),
  n_layers_by_core %>% transmute(metric = paste0("n_layers_", core_clean), value = n_layers, note = "Layer rows per core"),
  tibble(metric = "pseudocount", value = pseudocount, note = "Pseudocount added before log10 transform"),
  key_pooled %>%
    transmute(
      metric = paste0("pooled_estimate_", ratio_name),
      value = estimate,
      note = paste0("Pooled model standardized coefficient; direction=", estimate_direction)
    ),
  key_pooled %>%
    transmute(
      metric = paste0("pooled_p_", ratio_name),
      value = p_value,
      note = "Pooled model p-value"
    )
)

methods_lines <- c(
  "Supplementary analysis: XRF ratios vs library_concentration",
  "",
  "Data inputs:",
  "- data/metadata_v5.tsv",
  "- data/combined_xrf_not_normalized.tsv",
  "",
  "Layer construction and matching:",
  "- Core names were normalized by stripping _R1/_R2 suffixes to define core_clean.",
  "- Technical replicates were collapsed to one layer row per core_clean + depth_in_core_cm + y_bp.",
  "- Main response: library_concentration_mean and log10(library_concentration_mean + pseudocount).",
  "- XRF rows were aggregated to layer depth bins defined by adjacent midpoint boundaries; per-layer ratio value was the interval median.",
  "- Interval coverage diagnostics retained n_xrf_rows_interval and n_valid_<ratio>.",
  "",
  "XRF ratio framework:",
  "- Ratios were derived from non-normalized XRF intensities and interpreted only relatively.",
  "- Required ratios: ba_ti, p_ti, br_ti, si_ti, ca_ti, al_si, ti_al, fe_al, rb_sr, mn_fe.",
  "- Optional sensitivity ratios available: zr_al, rb_al.",
  "",
  "Core-specific analyses:",
  "- Raw Spearman and detrended Spearman (residualized against age or depth trend) were computed per ratio.",
  "- Per-ratio depth/age-adjusted models were fit with spline trend terms and standardized ratio effects.",
  "- ST13 and GeoB25202 attempted autocorrelation-aware GLS (corCAR1); fallback model was linear regression when GLS was not feasible.",
  "",
  "Cross-core pooled analyses:",
  "- Residual correlations adjusted for core and smooth age/depth trend.",
  "- Pooled per-ratio regression: response ~ ratio + core + spline(trend).",
  "- Interaction models (ratio x core) were fit for rb_sr and ca_ti to quantify heterogeneity.",
  "",
  "Interpretation guardrails:",
  "- Associations are observational and cannot establish mechanism.",
  "- Library concentration reflects DNA recovery and may also capture preservation and extractability effects.",
  "- Age/depth structure and core identity are major confounders.",
  "- Within-core results are primary for cautious mechanistic interpretation.",
  "- Pooled results are used to identify broad directional patterns and heterogeneity, not to erase core-specific differences."
)

notes_lines <- c(
  "Supplementary notes: library concentration and unified ratio-only XRF analysis",
  "",
  "Conservative interpretation:",
  "- Ratios were derived from non-normalized XRF and interpreted relatively only.",
  "- No direct interpretation of raw XRF elemental counts was used for scientific conclusions.",
  "- Library concentration is a DNA recovery proxy that can integrate preservation and extraction-related variation.",
  "- Age/depth gradients and core identity are expected confounders; all inferential summaries account for this structure.",
  "- Within-core estimates are prioritized for mechanistic caution.",
  "- Cross-core pooling is used only for broad directional tendencies and heterogeneity assessment.",
  "",
  "Model limitations:",
  "- Sample sizes per core can limit interaction precision.",
  "- Spline trend adjustment is flexible but does not eliminate all potential temporal confounding.",
  "- Ratio collinearity can affect single-ratio model comparability across cores."
)

# -----------------------------------------------------------------------------
# 10) Figures
# -----------------------------------------------------------------------------
ratio_order <- required_ratios

fig_core_forest <- core_models %>%
  filter(is.finite(estimate)) %>%
  mutate(
    ratio_name = factor(ratio_name, levels = rev(ratio_order)),
    core_clean = factor(core_clean, levels = params$core_targets)
  ) %>%
  ggplot(aes(x = estimate, y = ratio_name, color = ratio_family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.2, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 2.0) +
  facet_wrap(~core_clean, scales = "free_y") +
  labs(
    title = "Per-core depth/age-adjusted ratio effects on log10 library concentration",
    subtitle = "Standardized coefficients from ratio-specific models",
    x = "Coefficient (ratio_z)", y = "XRF ratio", color = "Ratio family"
  ) +
  theme_bw()
save_plot_both(fig_core_forest, "library_concentration_xrf_supplementary_per_core_forest", width = 12, height = 8)

fig_pooled <- pooled_models %>%
  filter(is.finite(estimate)) %>%
  mutate(ratio_name = factor(ratio_name, levels = rev(ratio_order))) %>%
  ggplot(aes(x = estimate, y = ratio_name, color = ratio_family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.2, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 2.2) +
  labs(
    title = "Cross-core pooled ratio effects",
    subtitle = "Model: log10(library concentration) ~ ratio + core + spline(trend)",
    x = "Pooled standardized coefficient", y = "XRF ratio", color = "Ratio family"
  ) +
  theme_bw()
save_plot_both(fig_pooled, "library_concentration_xrf_supplementary_pooled_effects", width = 9, height = 6)

if (nrow(interaction_detail) > 0) {
  fig_heterogeneity <- interaction_detail %>%
    mutate(core_clean = factor(core_clean, levels = params$core_targets)) %>%
    ggplot(aes(x = estimate, y = core_clean, color = ratio_name)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.2, alpha = 0.8) +
    geom_point(size = 2.2) +
    facet_wrap(~ratio_name, scales = "free_x") +
    labs(
      title = "Key-ratio heterogeneity across cores",
      subtitle = "Core-specific slopes from pooled ratio × core interaction models",
      x = "Core-specific standardized slope", y = "Core", color = "Ratio"
    ) +
    theme_bw()
  save_plot_both(fig_heterogeneity, "library_concentration_xrf_supplementary_key_ratio_heterogeneity", width = 10, height = 6)
}

fig_geob <- core_models %>%
  filter(core_clean == "GeoB25202", is.finite(estimate)) %>%
  mutate(ratio_name = factor(ratio_name, levels = rev(ratio_order))) %>%
  ggplot(aes(x = estimate, y = ratio_name, color = ratio_family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.2, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 2.4) +
  labs(
    title = "GeoB25202 focused ratio structure",
    subtitle = "Unified ratio-only XRF framework reveals multi-ratio gradient patterns",
    x = "Standardized coefficient", y = "XRF ratio", color = "Ratio family"
  ) +
  theme_bw()
save_plot_both(fig_geob, "library_concentration_xrf_supplementary_geob25202_focus", width = 8, height = 6)

fig_st13 <- core_models %>%
  filter(core_clean == "ST13", is.finite(estimate)) %>%
  mutate(ratio_name = factor(ratio_name, levels = rev(ratio_order))) %>%
  ggplot(aes(x = estimate, y = ratio_name, color = ratio_family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.2, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 2.4) +
  labs(
    title = "ST13 focused ratio structure",
    subtitle = "Strong terrigenous/mineralogical signal with accompanying productivity/carbonate terms",
    x = "Standardized coefficient", y = "XRF ratio", color = "Ratio family"
  ) +
  theme_bw()
save_plot_both(fig_st13, "library_concentration_xrf_supplementary_st13_focus", width = 8, height = 6)

fig_family <- core_models %>%
  filter(is.finite(estimate)) %>%
  mutate(direction = ifelse(estimate > 0, "positive", "negative")) %>%
  count(ratio_family, direction, name = "n_coeff") %>%
  ggplot(aes(x = ratio_family, y = n_coeff, fill = direction)) +
  geom_col(position = "stack") +
  labs(
    title = "Ratio-family direction synthesis",
    subtitle = "Count of per-core coefficient directions",
    x = "Ratio family", y = "Count of coefficients", fill = "Direction"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot_both(fig_family, "library_concentration_xrf_supplementary_ratio_family_directions", width = 8, height = 5)

# -----------------------------------------------------------------------------
# 11) Write outputs
# -----------------------------------------------------------------------------
write_tsv(layer_table, file.path(out_dir, "library_concentration_xrf_supplementary_layer_table.tsv"))
write_tsv(qc_summary, file.path(out_dir, "library_concentration_xrf_supplementary_qc_summary.tsv"))
write_tsv(core_associations, file.path(out_dir, "library_concentration_xrf_supplementary_core_associations.tsv"))
write_tsv(core_models, file.path(out_dir, "library_concentration_xrf_supplementary_core_models.tsv"))
write_tsv(pooled_models, file.path(out_dir, "library_concentration_xrf_supplementary_pooled_models.tsv"))
write_tsv(ratio_family_summary, file.path(out_dir, "library_concentration_xrf_supplementary_ratio_family_summary.tsv"))
write_tsv(interaction_detail, file.path(out_dir, "library_concentration_xrf_supplementary_interaction_summary.tsv"))
write_tsv(key_numbers, file.path(out_dir, "library_concentration_xrf_supplementary_key_numbers.tsv"))

writeLines(notes_lines, con = file.path(out_dir, "library_concentration_xrf_supplementary_notes.txt"))
writeLines(methods_lines, con = file.path(out_dir, "library_concentration_xrf_supplementary_methods.txt"))

# Keep pooled residual diagnostics in pooled models table directory context via message.
write_tsv(
  pooled_residual_correlations,
  file.path(out_dir, "library_concentration_xrf_supplementary_pooled_residual_correlations.tsv")
)

message("Supplementary XRF-library concentration analysis complete.")
message("Outputs: ", out_dir)
message("Plots: ", plot_dir)
