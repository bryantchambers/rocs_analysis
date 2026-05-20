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

# ------------------------------------------------------------------
# Paths, outputs, and parameters
# ------------------------------------------------------------------
metadata_path <- file.path("data", "metadata_v5.tsv")
xrf_path <- file.path("data", "combined_xrf_not_normalized.tsv")

out_dir <- file.path("results", "common", "proxies", "final")
plot_dir <- file.path("plots", "proxies", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(metadata_path, xrf_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files (run from repo root): ", paste(missing_files, collapse = ", "))
}

params <- list(
  min_n_correlation = 8,
  min_n_model = 20,
  trend_spline_df = 4,
  model_p_threshold = 0.10,
  pseudocount_rule = "min_positive_library_concentration_mean / 2",
  layer_key = "core_clean + depth_in_core_cm + y_bp",
  xrf_aggregation = "median within depth-interval bounds",
  response_primary = "log10(library_concentration_mean + pseudocount)"
)

ratio_definitions <- tibble::tribble(
  ~ratio_name, ~numerator, ~denominator, ~ratio_family, ~priority_group,
  "ba_ti", "ba", "ti", "productivity_biogenic", "required",
  "p_ti",  "p",  "ti", "productivity_biogenic", "required",
  "br_ti", "br", "ti", "productivity_biogenic", "required",
  "si_ti", "si", "ti", "productivity_biogenic", "required",
  "ca_ti", "ca", "ti", "productivity_biogenic", "required",
  "al_si", "al", "si", "terrigenous_mineralogical", "required",
  "ti_al", "ti", "al", "terrigenous_mineralogical", "required",
  "fe_al", "fe", "al", "terrigenous_mineralogical", "required",
  "rb_sr", "rb", "sr", "terrigenous_mineralogical", "required",
  "mn_fe", "mn", "fe", "redox", "required",
  "zr_al", "zr", "al", "terrigenous_mineralogical", "optional",
  "rb_al", "rb", "al", "terrigenous_mineralogical", "optional"
)

core_model_targets <- c("ST13", "ST5", "ST8", "GeoB25202")

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
normalize_core <- function(x) str_replace(x, "_R[12]$", "")

safe_min_positive <- function(x, default = 1e-6) {
  v <- x[is.finite(x) & x > 0]
  if (length(v) == 0) return(default)
  min(v)
}

safe_cor <- function(x, y, min_n = 8, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < min_n || dplyr::n_distinct(x2) < 3 || dplyr::n_distinct(y2) < 3) {
    return(list(n = length(x2), estimate = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = method, exact = FALSE))
  list(n = length(x2), estimate = unname(ct$estimate), p = ct$p.value)
}

choose_trend_col <- function(df) {
  if ("y_bp" %in% names(df)) {
    ok <- is.finite(df$y_bp)
    if (sum(ok) >= params$min_n_correlation && dplyr::n_distinct(df$y_bp[ok]) >= 6) return("y_bp")
  }
  if ("depth_in_core_cm" %in% names(df)) {
    ok <- is.finite(df$depth_in_core_cm)
    if (sum(ok) >= params$min_n_correlation && dplyr::n_distinct(df$depth_in_core_cm[ok]) >= 6) return("depth_in_core_cm")
  }
  NA_character_
}

safe_ns_residuals <- function(v, trend, target_df = 4) {
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v) & is.finite(trend)
  if (sum(ok) < params$min_n_correlation) return(out)

  k <- dplyr::n_distinct(trend[ok])
  df_use <- max(1, min(target_df, k - 1))
  fit <- if (k >= 4) {
    tryCatch(lm(v[ok] ~ splines::ns(trend[ok], df = df_use)), error = function(e) NULL)
  } else {
    tryCatch(lm(v[ok] ~ trend[ok]), error = function(e) NULL)
  }
  if (is.null(fit)) return(out)
  out[ok] <- residuals(fit)
  out
}

layer_intervals <- function(depths) {
  d <- sort(unique(depths[is.finite(depths)]))
  if (length(d) == 0) {
    return(tibble(depth_in_core_cm = numeric(), lower = numeric(), upper = numeric()))
  }
  if (length(d) == 1) {
    return(tibble(depth_in_core_cm = d, lower = -Inf, upper = Inf))
  }
  mids <- (d[-1] + d[-length(d)]) / 2
  tibble(depth_in_core_cm = d, lower = c(-Inf, mids), upper = c(mids, Inf))
}

derive_ratio <- function(df, num_col, den_col, out_col) {
  if (!(num_col %in% names(df)) || !(den_col %in% names(df))) {
    df[[out_col]] <- NA_real_
    return(df)
  }
  num <- suppressWarnings(as.numeric(df[[num_col]]))
  den <- suppressWarnings(as.numeric(df[[den_col]]))
  df[[out_col]] <- ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
  df
}

aggregate_xrf_to_layers <- function(layer_df, xrf_core_df, ratio_names) {
  ints <- layer_intervals(layer_df$depth_in_core_cm)
  if (nrow(ints) == 0) {
    out <- layer_df %>% select(depth_in_core_cm)
    out$n_xrf_rows_interval <- integer(nrow(out))
    for (rn in ratio_names) {
      out[[rn]] <- NA_real_
      out[[paste0("n_valid_", rn)]] <- NA_integer_
    }
    return(out)
  }

  out <- ints %>% select(depth_in_core_cm)
  out$n_xrf_rows_interval <- NA_integer_
  for (rn in ratio_names) {
    out[[rn]] <- NA_real_
    out[[paste0("n_valid_", rn)]] <- NA_integer_
  }

  for (i in seq_len(nrow(ints))) {
    lo <- ints$lower[i]
    hi <- ints$upper[i]
    idx <- which(
      is.finite(xrf_core_df$depth_in_core_cm) &
        xrf_core_df$depth_in_core_cm > lo &
        xrf_core_df$depth_in_core_cm <= hi
    )
    out$n_xrf_rows_interval[i] <- length(idx)
    if (length(idx) == 0) next

    for (rn in ratio_names) {
      vv <- xrf_core_df[[rn]][idx]
      v_ok <- vv[is.finite(vv)]
      out[[paste0("n_valid_", rn)]][i] <- length(v_ok)
      if (length(v_ok) > 0) {
        out[[rn]][i] <- median(v_ok)
      }
    }
  }

  out
}

compute_core_ratio_associations <- function(df_core, ratio_info) {
  trend_col <- choose_trend_col(df_core)

  map_dfr(seq_len(nrow(ratio_info)), function(i) {
    rn <- ratio_info$ratio_name[i]
    fam <- ratio_info$ratio_family[i]
    pri <- ratio_info$priority_group[i]

    x <- df_core[[rn]]
    y <- df_core$library_concentration_log10
    raw <- safe_cor(x, y, min_n = params$min_n_correlation)

    det <- list(n = NA_integer_, estimate = NA_real_, p = NA_real_)
    if (!is.na(trend_col)) {
      xr <- safe_ns_residuals(x, df_core[[trend_col]], target_df = params$trend_spline_df)
      yr <- safe_ns_residuals(y, df_core[[trend_col]], target_df = params$trend_spline_df)
      det <- safe_cor(xr, yr, min_n = params$min_n_correlation)
    }

    tibble(
      core_clean = unique(df_core$core_clean)[1],
      ratio_name = rn,
      ratio_family = fam,
      priority_group = pri,
      trend_used = trend_col,
      n_raw = raw$n,
      spearman_raw = raw$estimate,
      p_raw = raw$p,
      n_detrended = det$n,
      spearman_detrended = det$estimate,
      p_detrended = det$p
    )
  })
}

fit_ratio_model <- function(df, ratio_name, core_adjusted = FALSE, try_gls = FALSE) {
  trend_col <- choose_trend_col(df)
  if (is.na(trend_col)) {
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

  use_cols <- c("library_concentration_log10", trend_col, ratio_name)
  if (core_adjusted) use_cols <- c(use_cols, "core_clean")

  dat <- df %>%
    select(all_of(use_cols)) %>%
    filter(is.finite(library_concentration_log10), is.finite(.data[[trend_col]]), is.finite(.data[[ratio_name]]))

  if (nrow(dat) < params$min_n_model || dplyr::n_distinct(dat[[ratio_name]]) < 4) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_col,
      model_type = "not_fit",
      n_model = nrow(dat),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "too few rows or ratio variation"
    ))
  }

  dat <- dat %>%
    mutate(
      ratio_z = as.numeric(scale(.data[[ratio_name]]))
    )

  if (!all(is.finite(dat$ratio_z))) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_col,
      model_type = "not_fit",
      n_model = nrow(dat),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "ratio standardization failed"
    ))
  }

  k <- dplyr::n_distinct(dat[[trend_col]])
  df_ns <- max(1, min(params$trend_spline_df, k - 1))

  rhs <- paste0("splines::ns(", trend_col, ", df = ", df_ns, ") + ratio_z")
  if (core_adjusted) rhs <- paste(rhs, "+ core_clean")
  fml <- as.formula(paste("library_concentration_log10 ~", rhs))

  if (try_gls && requireNamespace("nlme", quietly = TRUE)) {
    dat <- dat %>% arrange(.data[[trend_col]])
    gls_fit <- tryCatch(
      nlme::gls(
        model = fml,
        data = dat,
        method = "REML",
        correlation = nlme::corCAR1(form = stats::as.formula(paste0("~", trend_col)))
      ),
      error = function(e) NULL
    )
    if (!is.null(gls_fit)) {
      tt <- summary(gls_fit)$tTable
      if ("ratio_z" %in% rownames(tt)) {
        ci <- tryCatch(nlme::intervals(gls_fit, which = "coef")$coef, error = function(e) NULL)
        return(tibble(
          ratio_name = ratio_name,
          trend_used = trend_col,
          model_type = "gls_corCAR1",
          n_model = nrow(dat),
          estimate = tt["ratio_z", "Value"],
          std_error = tt["ratio_z", "Std.Error"],
          p_value = tt["ratio_z", "p-value"],
          ci_low = if (!is.null(ci)) ci["ratio_z", "lower"] else NA_real_,
          ci_high = if (!is.null(ci)) ci["ratio_z", "upper"] else NA_real_,
          note = paste0("depth/age trend ns df=", df_ns)
        ))
      }
    }
  }

  lm_fit <- tryCatch(lm(fml, data = dat), error = function(e) NULL)
  if (is.null(lm_fit)) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_col,
      model_type = "not_fit",
      n_model = nrow(dat),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      note = "lm fit error"
    ))
  }

  sm <- summary(lm_fit)$coefficients
  if (!("ratio_z" %in% rownames(sm))) {
    return(tibble(
      ratio_name = ratio_name,
      trend_used = trend_col,
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
    trend_used = trend_col,
    model_type = "lm",
    n_model = nrow(dat),
    estimate = sm["ratio_z", "Estimate"],
    std_error = sm["ratio_z", "Std. Error"],
    p_value = sm["ratio_z", "Pr(>|t|)"],
    ci_low = if (!is.null(ci)) ci["ratio_z", 1] else NA_real_,
    ci_high = if (!is.null(ci)) ci["ratio_z", 2] else NA_real_,
    note = paste0("depth/age trend ns df=", df_ns)
  )
}

# ------------------------------------------------------------------
# 1) Load and collapse metadata to sediment-layer rows
# ------------------------------------------------------------------
metadata_raw <- read_tsv(metadata_path, show_col_types = FALSE)
xrf_raw <- read_tsv(xrf_path, show_col_types = FALSE)

metadata_layers <- metadata_raw %>%
  mutate(core_clean = normalize_core(core)) %>%
  filter(
    core_clean %in% unique(xrf_raw$core),
    is.finite(depth_in_core_cm),
    is.finite(y_bp)
  ) %>%
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

# ------------------------------------------------------------------
# 2) Derive ratio-only proxies from non-normalized XRF table
# ------------------------------------------------------------------
xrf <- xrf_raw %>%
  mutate(core_clean = as.character(core))

for (i in seq_len(nrow(ratio_definitions))) {
  xrf <- derive_ratio(
    xrf,
    num_col = ratio_definitions$numerator[i],
    den_col = ratio_definitions$denominator[i],
    out_col = ratio_definitions$ratio_name[i]
  )
}

ratio_names <- ratio_definitions$ratio_name

# ------------------------------------------------------------------
# 3) Aggregate dense XRF to library-layer intervals
# ------------------------------------------------------------------
agg_layers <- map_dfr(unique(metadata_layers$core_clean), function(cc) {
  layer_core <- metadata_layers %>% filter(core_clean == cc)
  xrf_core <- xrf %>% filter(core_clean == cc, is.finite(depth_in_core_cm))

  if (nrow(layer_core) == 0) return(tibble())

  agg <- aggregate_xrf_to_layers(layer_core, xrf_core, ratio_names = ratio_names)
  layer_core %>%
    left_join(agg, by = "depth_in_core_cm")
})

analysis_table <- agg_layers %>%
  mutate(core_clean = as.character(core_clean))

# ------------------------------------------------------------------
# 4) QC summary
# ------------------------------------------------------------------
qc_summary <- map_dfr(unique(analysis_table$core_clean), function(cc) {
  core_xrf <- xrf %>% filter(core_clean == cc)
  core_layers <- analysis_table %>% filter(core_clean == cc)

  map_dfr(seq_len(nrow(ratio_definitions)), function(i) {
    rn <- ratio_definitions$ratio_name[i]
    num <- ratio_definitions$numerator[i]
    den <- ratio_definitions$denominator[i]

    den_vals <- if (den %in% names(core_xrf)) suppressWarnings(as.numeric(core_xrf[[den]])) else rep(NA_real_, nrow(core_xrf))
    num_vals <- if (num %in% names(core_xrf)) suppressWarnings(as.numeric(core_xrf[[num]])) else rep(NA_real_, nrow(core_xrf))

    tibble(
      core_clean = cc,
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[i],
      priority_group = ratio_definitions$priority_group[i],
      numerator = num,
      denominator = den,
      ratio_derivable_in_source = any(is.finite(core_xrf[[rn]])),
      n_xrf_rows_core = nrow(core_xrf),
      n_layers_core = nrow(core_layers),
      n_layers_with_response = sum(is.finite(core_layers$library_concentration_log10)),
      n_layers_with_ratio = sum(is.finite(core_layers[[rn]])),
      pct_layers_with_ratio = 100 * mean(is.finite(core_layers[[rn]])),
      median_n_xrf_rows_interval = median(core_layers$n_xrf_rows_interval, na.rm = TRUE),
      median_n_valid_ratio_per_interval = median(core_layers[[paste0("n_valid_", rn)]], na.rm = TRUE),
      denominator_nonpositive_or_missing_rows = sum(!is.finite(den_vals) | den_vals <= 0),
      numerator_missing_rows = sum(!is.finite(num_vals))
    )
  })
}) %>%
  mutate(
    median_n_xrf_rows_interval = ifelse(is.nan(median_n_xrf_rows_interval), NA_real_, median_n_xrf_rows_interval),
    median_n_valid_ratio_per_interval = ifelse(is.nan(median_n_valid_ratio_per_interval), NA_real_, median_n_valid_ratio_per_interval)
  ) %>%
  arrange(core_clean, ratio_family, ratio_name)

# ------------------------------------------------------------------
# 5) Core-specific exploratory correlations (raw + detrended)
# ------------------------------------------------------------------
assoc_core <- map_dfr(unique(analysis_table$core_clean), function(cc) {
  compute_core_ratio_associations(
    df_core = analysis_table %>% filter(core_clean == cc),
    ratio_info = ratio_definitions
  )
})

assoc_core <- assoc_core %>%
  group_by(core_clean) %>%
  mutate(
    q_raw_bh_core = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh_core = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  ungroup() %>%
  group_by(core_clean, ratio_family) %>%
  mutate(
    q_raw_bh_core_family = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh_core_family = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  ungroup()

# ------------------------------------------------------------------
# 6) Pooled core-adjusted perspective (secondary)
# ------------------------------------------------------------------
pooled_models <- map_dfr(ratio_names, function(rn) {
  fit_ratio_model(
    df = analysis_table,
    ratio_name = rn,
    core_adjusted = TRUE,
    try_gls = FALSE
  ) %>%
    mutate(core_clean = "ALL_CORES_ADJUSTED", analysis_type = "pooled_depth_core_adjusted_lm")
})

pooled_residual_corr <- map_dfr(ratio_names, function(rn) {
  dat <- analysis_table %>%
    select(core_clean, y_bp, depth_in_core_cm, library_concentration_log10, all_of(rn)) %>%
    filter(is.finite(library_concentration_log10), is.finite(.data[[rn]]))

  trend_col <- choose_trend_col(dat)
  if (nrow(dat) < params$min_n_correlation || is.na(trend_col)) {
    return(tibble(
      core_clean = "ALL_CORES_ADJUSTED",
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
      priority_group = ratio_definitions$priority_group[match(rn, ratio_definitions$ratio_name)],
      trend_used = trend_col,
      n_raw = nrow(dat),
      spearman_raw = NA_real_,
      p_raw = NA_real_,
      n_detrended = nrow(dat),
      spearman_detrended = NA_real_,
      p_detrended = NA_real_,
      q_raw_bh_core = NA_real_,
      q_detrended_bh_core = NA_real_,
      q_raw_bh_core_family = NA_real_,
      q_detrended_bh_core_family = NA_real_,
      analysis_type = "pooled_core_plus_trend_residual_spearman"
    ))
  }

  k <- dplyr::n_distinct(dat[[trend_col]])
  df_ns <- max(1, min(params$trend_spline_df, k - 1))
  fm_y <- as.formula(paste0("library_concentration_log10 ~ splines::ns(", trend_col, ", df=", df_ns, ") + core_clean"))
  fm_x <- as.formula(paste0(rn, " ~ splines::ns(", trend_col, ", df=", df_ns, ") + core_clean"))

  fit_y <- tryCatch(lm(fm_y, data = dat), error = function(e) NULL)
  fit_x <- tryCatch(lm(fm_x, data = dat), error = function(e) NULL)

  if (is.null(fit_y) || is.null(fit_x)) {
    return(tibble(
      core_clean = "ALL_CORES_ADJUSTED",
      ratio_name = rn,
      ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
      priority_group = ratio_definitions$priority_group[match(rn, ratio_definitions$ratio_name)],
      trend_used = trend_col,
      n_raw = nrow(dat),
      spearman_raw = NA_real_,
      p_raw = NA_real_,
      n_detrended = nrow(dat),
      spearman_detrended = NA_real_,
      p_detrended = NA_real_,
      q_raw_bh_core = NA_real_,
      q_detrended_bh_core = NA_real_,
      q_raw_bh_core_family = NA_real_,
      q_detrended_bh_core_family = NA_real_,
      analysis_type = "pooled_core_plus_trend_residual_spearman"
    ))
  }

  y_res <- residuals(fit_y)
  x_res <- residuals(fit_x)
  det <- safe_cor(x_res, y_res, min_n = params$min_n_correlation)
  raw <- safe_cor(dat[[rn]], dat$library_concentration_log10, min_n = params$min_n_correlation)

  tibble(
    core_clean = "ALL_CORES_ADJUSTED",
    ratio_name = rn,
    ratio_family = ratio_definitions$ratio_family[match(rn, ratio_definitions$ratio_name)],
    priority_group = ratio_definitions$priority_group[match(rn, ratio_definitions$ratio_name)],
    trend_used = trend_col,
    n_raw = raw$n,
    spearman_raw = raw$estimate,
    p_raw = raw$p,
    n_detrended = det$n,
    spearman_detrended = det$estimate,
    p_detrended = det$p,
    q_raw_bh_core = NA_real_,
    q_detrended_bh_core = NA_real_,
    q_raw_bh_core_family = NA_real_,
    q_detrended_bh_core_family = NA_real_,
    analysis_type = "pooled_core_plus_trend_residual_spearman"
  )
}) %>%
  mutate(
    q_raw_bh_core = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh_core = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  group_by(ratio_family) %>%
  mutate(
    q_raw_bh_core_family = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh_core_family = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  ungroup()

# ------------------------------------------------------------------
# 7) Core-specific depth-adjusted model summaries
# ------------------------------------------------------------------
core_model_ratios <- ratio_definitions %>%
  filter(priority_group == "required")

core_models <- map_dfr(core_model_targets, function(cc) {
  dat <- analysis_table %>% filter(core_clean == cc)
  robust_gls <- cc %in% c("ST13", "GeoB25202")

  map_dfr(core_model_ratios$ratio_name, function(rn) {
    fit_ratio_model(
      df = dat,
      ratio_name = rn,
      core_adjusted = FALSE,
      try_gls = robust_gls
    ) %>%
      mutate(
        core_clean = cc,
        ratio_family = core_model_ratios$ratio_family[match(rn, core_model_ratios$ratio_name)],
        priority_group = core_model_ratios$priority_group[match(rn, core_model_ratios$ratio_name)],
        estimate_sign = case_when(
          is.na(estimate) ~ NA_character_,
          estimate > 0 ~ "positive",
          estimate < 0 ~ "negative",
          TRUE ~ "zero"
        )
      )
  })
}) %>%
  select(
    core_clean, ratio_name, ratio_family, priority_group,
    model_type, trend_used, n_model,
    estimate, estimate_sign, std_error, p_value, ci_low, ci_high, note
  ) %>%
  group_by(core_clean) %>%
  mutate(q_model_bh_core = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH"))) %>%
  ungroup() %>%
  group_by(core_clean, ratio_family) %>%
  mutate(q_model_bh_core_family = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH"))) %>%
  ungroup()

# ------------------------------------------------------------------
# 8) Compile unified association output
# ------------------------------------------------------------------
assoc_unified <- bind_rows(
  assoc_core %>% mutate(analysis_type = "within_core_spearman"),
  pooled_residual_corr
) %>%
  arrange(core_clean, ratio_family, ratio_name)

# ------------------------------------------------------------------
# 9) Prior-comparison logic
# ------------------------------------------------------------------
get_core_model_row <- function(core_id, ratio_id) {
  core_models %>%
    filter(core_clean == core_id, ratio_name == ratio_id) %>%
    slice_head(n = 1)
}

get_core_assoc_row <- function(core_id, ratio_id) {
  assoc_core %>%
    filter(core_clean == core_id, ratio_name == ratio_id) %>%
    slice_head(n = 1)
}

st13_rb <- get_core_model_row("ST13", "rb_sr")
st13_rb_assoc <- get_core_assoc_row("ST13", "rb_sr")

st13_rb_call <- if (nrow(st13_rb) == 1 && is.finite(st13_rb$estimate) && is.finite(st13_rb$p_value)) {
  if (st13_rb$estimate < 0 && st13_rb$p_value <= 0.05) "holds_and_stronger"
  else if (st13_rb$estimate < 0 && st13_rb$p_value <= 0.10) "holds"
  else if (st13_rb$estimate < 0) "holds_but_weaker"
  else "changed"
} else {
  "inconclusive"
}

geob_prod_ratios <- c("ca_ti", "ba_ti", "br_ti", "p_ti", "si_ti")
geob_prod <- core_models %>%
  filter(core_clean == "GeoB25202", ratio_name %in% geob_prod_ratios)

geob_prod_sig <- geob_prod %>% filter(is.finite(p_value), p_value <= 0.10)
geob_prod_call <- if (nrow(geob_prod_sig) >= 2) {
  "stronger_available_signal"
} else if (nrow(geob_prod_sig) == 1) {
  "some_new_signal"
} else {
  "still_limited"
}

st5_terr <- core_models %>% filter(core_clean == "ST5", ratio_family == "terrigenous_mineralogical")
st8_terr <- core_models %>% filter(core_clean == "ST8", ratio_family == "terrigenous_mineralogical")
st5_terr_support <- any(is.finite(st5_terr$p_value) & st5_terr$p_value <= 0.10 & st5_terr$estimate < 0)
st8_terr_support <- any(is.finite(st8_terr$p_value) & st8_terr$p_value <= 0.10 & st8_terr$estimate < 0)
st58_call <- if (!st5_terr_support && !st8_terr_support) "holds_largely_unsupported" else "mixed_or_changed"

terr_sig <- core_models %>%
  filter(ratio_family == "terrigenous_mineralogical") %>%
  summarise(n_sig = sum(is.finite(p_value) & p_value <= 0.10), .groups = "drop") %>%
  pull(n_sig)
prod_sig <- core_models %>%
  filter(ratio_family == "productivity_biogenic") %>%
  summarise(n_sig = sum(is.finite(p_value) & p_value <= 0.10), .groups = "drop") %>%
  pull(n_sig)

dominance_call <- if (is.na(terr_sig) || is.na(prod_sig)) {
  "inconclusive"
} else if (terr_sig > prod_sig) {
  "still_supported"
} else if (terr_sig == prod_sig) {
  "mixed"
} else {
  "weakened"
}

prior_comparison <- bind_rows(
  tibble(
    question_id = 1,
    prior_question = "Does the negative Rb/Sr signal still hold in ST13?",
    previous_conclusion = "Yes, negative and robust in prior analyses.",
    new_result_summary = paste0(
      "ST13 rb_sr model estimate=", formatC(st13_rb$estimate, format = "f", digits = 3),
      ", p=", formatC(st13_rb$p_value, format = "e", digits = 2),
      "; detrended rho=", formatC(st13_rb_assoc$spearman_detrended, format = "f", digits = 3),
      ", p=", formatC(st13_rb_assoc$p_detrended, format = "e", digits = 2)
    ),
    comparison_call = st13_rb_call,
    evidence_basis = "ST13 depth/age-adjusted core-specific model + detrended Spearman"
  ),
  tibble(
    question_id = 2,
    prior_question = "Does GeoB25202 now show productivity/carbonate signals from unified ratios?",
    previous_conclusion = "Prior GeoB25202 interpretation was more limited/heterogeneous by proxy source.",
    new_result_summary = paste0(
      "GeoB25202 productivity ratios with p<=0.10: ",
      ifelse(nrow(geob_prod_sig) == 0, "none", paste(geob_prod_sig$ratio_name, collapse = ", "))
    ),
    comparison_call = geob_prod_call,
    evidence_basis = "GeoB25202 core-specific models for ca_ti, ba_ti, br_ti, p_ti, si_ti"
  ),
  tibble(
    question_id = 3,
    prior_question = "Do ST5/ST8 remain largely unsupported for terrigenous signal?",
    previous_conclusion = "Prior ST5/ST8 support was weak/inconsistent.",
    new_result_summary = paste0(
      "ST5 terrigenous support=", st5_terr_support,
      "; ST8 terrigenous support=", st8_terr_support
    ),
    comparison_call = st58_call,
    evidence_basis = "ST5/ST8 core-specific terrigenous model terms"
  ),
  tibble(
    question_id = 4,
    prior_question = "Is terrigenous/mineralogical dominance over productivity/carbonate still supported?",
    previous_conclusion = "Prior synthesis favored terrigenous/mineralogical context.",
    new_result_summary = paste0(
      "Count of p<=0.10 model terms: terrigenous=", terr_sig,
      "; productivity=", prod_sig
    ),
    comparison_call = dominance_call,
    evidence_basis = "Cross-core count of depth-adjusted model signals by ratio family"
  )
) %>%
  arrange(question_id)

# ------------------------------------------------------------------
# 10) Plots
# ------------------------------------------------------------------
cross_core_ratios <- c("rb_sr", "ca_ti", "ba_ti", "mn_fe")

plot_cross_core <- core_models %>%
  filter(ratio_name %in% cross_core_ratios, is.finite(estimate)) %>%
  mutate(
    ratio_name = factor(ratio_name, levels = cross_core_ratios),
    core_clean = factor(core_clean, levels = core_model_targets)
  )

if (nrow(plot_cross_core) > 0) {
  p1 <- ggplot(plot_cross_core, aes(x = estimate, y = core_clean, color = ratio_name)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(position = position_dodge(width = 0.5), size = 2.2) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      position = position_dodge(width = 0.5),
      width = 0.15,
      na.rm = TRUE
    ) +
    labs(
      title = "Cross-core depth-adjusted ratio effects on log library concentration",
      subtitle = "Ratios from non-normalized XRF source; observational associations only",
      x = "Standardized model coefficient (ratio_z)",
      y = "Core",
      color = "Ratio"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(plot_dir, "library_concentration_xrf_ratio_cross_core_effects.png"),
    plot = p1,
    width = 10,
    height = 6,
    dpi = 300
  )
}

geob_focus_ratios <- c("rb_sr", "ca_ti", "ba_ti", "br_ti", "p_ti", "si_ti")
plot_geob <- analysis_table %>%
  filter(core_clean == "GeoB25202") %>%
  select(core_clean, depth_in_core_cm, y_bp, library_concentration_log10, all_of(intersect(geob_focus_ratios, names(.)))) %>%
  pivot_longer(
    cols = all_of(intersect(geob_focus_ratios, names(.))),
    names_to = "ratio_name",
    values_to = "ratio_value"
  ) %>%
  filter(is.finite(ratio_value), is.finite(library_concentration_log10))

if (nrow(plot_geob) > 0) {
  p2 <- ggplot(plot_geob, aes(x = ratio_value, y = library_concentration_log10)) +
    geom_point(size = 1.3, alpha = 0.7, color = "#1B5E20") +
    geom_smooth(method = "lm", se = FALSE, color = "#0D47A1", linewidth = 0.8) +
    facet_wrap(~ratio_name, scales = "free_x", ncol = 3) +
    labs(
      title = "GeoB25202: library concentration vs key unified XRF ratios",
      subtitle = "New non-normalized XRF table (ratio-only interpretation)",
      x = "Ratio value",
      y = "log10(library_concentration_mean + pseudocount)"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(plot_dir, "library_concentration_xrf_ratio_geob25202_focus.png"),
    plot = p2,
    width = 12,
    height = 8,
    dpi = 300
  )
}

st13_focus <- assoc_core %>%
  filter(core_clean == "ST13", ratio_name %in% ratio_definitions$ratio_name) %>%
  mutate(
    ratio_name = factor(ratio_name, levels = ratio_definitions$ratio_name),
    ratio_family = factor(ratio_family, levels = c("terrigenous_mineralogical", "productivity_biogenic", "redox"))
  )

if (nrow(st13_focus) > 0) {
  p3 <- ggplot(st13_focus, aes(x = spearman_detrended, y = ratio_name, color = ratio_family)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(size = 2.3) +
    labs(
      title = "ST13 detrended ratio associations by proxy family",
      subtitle = "Depth/age-adjusted Spearman rho",
      x = "Detrended Spearman rho",
      y = "Ratio",
      color = "Family"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(plot_dir, "library_concentration_xrf_ratio_st13_focus.png"),
    plot = p3,
    width = 9,
    height = 6,
    dpi = 300
  )
}

# ------------------------------------------------------------------
# 11) Notes and write outputs
# ------------------------------------------------------------------
notes_lines <- c(
  "Harmonized library_concentration vs unified XRF ratio re-analysis",
  "",
  "Interpretation guardrails:",
  "- Ratios were derived from a non-normalized XRF table and are interpreted only as relative proxies.",
  "- Raw elemental counts were not interpreted directly.",
  "- Associations are observational and not causal.",
  "- Depth/age structure and core identity remain major confounders.",
  "- Within-core interpretation is primary; pooled cross-core results are secondary and exploratory.",
  "",
  "Matching and aggregation rules:",
  "- Core names normalized by stripping _R1/_R2 to core_clean.",
  paste0("- Metadata technical replicates collapsed by ", params$layer_key, "."),
  "- Dense XRF aggregated to library layers by midpoint-defined depth intervals.",
  paste0("- Interval summary used: ", params$xrf_aggregation, "."),
  "- Number of XRF rows and valid-ratio rows per interval retained.",
  "",
  "Modeling notes:",
  paste0("- Main response: ", params$response_primary, "."),
  paste0("- Pseudocount rule: ", params$pseudocount_rule, "."),
  "- Core-specific depth/age-adjusted models fitted for required ratios.",
  "- For ST13 and GeoB25202, GLS with corCAR1 was attempted as an autocorrelation-robust fallback when feasible."
)

writeLines(notes_lines, con = file.path(out_dir, "library_concentration_xrf_ratio_notes.txt"))

write_tsv(qc_summary, file.path(out_dir, "library_concentration_xrf_ratio_qc_summary.tsv"))
write_tsv(assoc_unified, file.path(out_dir, "library_concentration_xrf_ratio_associations.tsv"))
write_tsv(core_models, file.path(out_dir, "library_concentration_xrf_ratio_core_models.tsv"))
write_tsv(prior_comparison, file.path(out_dir, "library_concentration_xrf_ratio_prior_comparison.tsv"))

message("Done. Wrote outputs to: ", out_dir)
message("Plots written to: ", plot_dir)
