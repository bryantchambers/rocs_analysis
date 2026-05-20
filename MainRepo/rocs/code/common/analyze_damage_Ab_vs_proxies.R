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
# Paths and output dirs
# ------------------------------------------------------------------
out_dir <- file.path("results", "common", "proxies", "final")
plot_dir <- file.path("plots", "proxies", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

damage_path <- file.path(
  "results", "microbial", "damage", "damage-classification-depositional", "sample-baselines.tsv"
)
classification_path_primary <- file.path(
  "results", "microbial", "damage", "damage-classification-depositional", "sample-classification-stats.tsv"
)
classification_path_fallback <- file.path(
  "results", "microbial", "damage", "damage-classification-depositional", "exploration", "sample-classification-stats.tsv"
)
metadata_path <- file.path("data", "metadata_v5.tsv")
xrf_path <- file.path("data", "combined_xrf_not_normalized.tsv")
foram_path <- file.path("data", "combined_foraminifera_geochem.tsv")
sst_path <- file.path("data", "combined_sst_proxies_separate_columns.tsv")
geob_path <- file.path("data", "geob25202_clean_proxies.tsv")

classification_path <- if (file.exists(classification_path_primary)) {
  classification_path_primary
} else {
  classification_path_fallback
}

required <- c(damage_path, classification_path, metadata_path, xrf_path, foram_path, sst_path, geob_path)
missing_files <- required[!file.exists(required)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------
min_rows_primary <- 100
sparse_depth_tolerance_cm <- 5
sparse_min_match_fraction <- 0.6
sparse_max_median_depth_diff_cm <- 5
min_n_correlation <- 8
min_n_model <- 25

normalize_core <- function(x) {
  str_replace(x, "_R[0-9]+$", "")
}

safe_cor <- function(x, y, method = "spearman", min_n = min_n_correlation) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < min_n || dplyr::n_distinct(x2) < 3 || dplyr::n_distinct(y2) < 3) {
    return(list(n = length(x2), estimate = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = method, exact = FALSE))
  list(n = length(x2), estimate = unname(ct$estimate), p = ct$p.value)
}

safe_ns_residuals <- function(v, trend, target_df = 4) {
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v) & is.finite(trend)
  if (sum(ok) < min_n_correlation) return(out)
  k <- length(unique(trend[ok]))
  if (k < 4) {
    fit <- tryCatch(lm(v[ok] ~ trend[ok]), error = function(e) NULL)
  } else {
    df_use <- max(1, min(target_df, k - 1))
    fit <- tryCatch(lm(v[ok] ~ splines::ns(trend[ok], df = df_use)), error = function(e) NULL)
  }
  if (is.null(fit)) return(out)
  out[ok] <- residuals(fit)
  out
}

choose_trend <- function(df) {
  if ("y_bp" %in% names(df)) {
    ok <- is.finite(df$y_bp)
    if (sum(ok) >= min_n_correlation && dplyr::n_distinct(df$y_bp[ok]) >= 4) return("y_bp")
  }
  if ("depth_in_core_cm" %in% names(df)) {
    ok <- is.finite(df$depth_in_core_cm)
    if (sum(ok) >= min_n_correlation && dplyr::n_distinct(df$depth_in_core_cm[ok]) >= 4) return("depth_in_core_cm")
  }
  NA_character_
}

ratio_or_na <- function(num, den) {
  ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
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

ratio_definitions <- tibble::tribble(
  ~ratio_name, ~numerator, ~denominator, ~ratio_family, ~priority_group,
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

layer_intervals <- function(depths) {
  d <- sort(unique(depths[is.finite(depths)]))
  if (length(d) == 0) {
    return(tibble(depth_in_core_cm = numeric(), lower = numeric(), upper = numeric()))
  }
  if (length(d) == 1) {
    return(tibble(depth_in_core_cm = d, lower = -Inf, upper = Inf))
  }
  mids <- (d[-1] + d[-length(d)]) / 2
  tibble(
    depth_in_core_cm = d,
    lower = c(-Inf, mids),
    upper = c(mids, Inf)
  )
}

aggregate_dense_to_layers <- function(layer_df, dense_df, proxy_cols) {
  intervals <- layer_intervals(layer_df$depth_in_core_cm)
  if (nrow(intervals) == 0 || nrow(dense_df) == 0) {
    out <- intervals %>% select(depth_in_core_cm)
    out$n_dense_rows <- integer(nrow(out))
    for (nm in proxy_cols) out[[nm]] <- NA_real_
    return(out)
  }

  out <- intervals %>% select(depth_in_core_cm)
  out$n_dense_rows <- NA_integer_
  for (nm in proxy_cols) out[[nm]] <- NA_real_

  for (i in seq_len(nrow(intervals))) {
    lo <- intervals$lower[i]
    hi <- intervals$upper[i]
    idx <- which(
      is.finite(dense_df$depth_in_core_cm) &
        dense_df$depth_in_core_cm > lo &
        dense_df$depth_in_core_cm <= hi
    )
    out$n_dense_rows[i] <- length(idx)
    if (length(idx) == 0) next
    for (nm in proxy_cols) {
      vv <- dense_df[[nm]][idx]
      vv <- vv[is.finite(vv)]
      if (length(vv) == 0) next
      out[[nm]][i] <- median(vv)
    }
  }
  out
}

nearest_one_row <- function(layer_row, proxy_df) {
  if (nrow(proxy_df) == 0) return(NULL)
  age_diff <- abs(proxy_df$y_bp - layer_row$y_bp)
  depth_diff <- abs(proxy_df$depth_in_core_cm - layer_row$depth_in_core_cm)

  age_diff[!is.finite(age_diff)] <- Inf
  depth_diff[!is.finite(depth_diff)] <- Inf
  idx <- which.min(age_diff)
  if (!is.finite(age_diff[idx])) {
    idx <- which.min(depth_diff)
  }
  if (!is.finite(depth_diff[idx]) && !is.finite(age_diff[idx])) return(NULL)

  proxy_df[idx, , drop = FALSE] %>%
    mutate(
      geob_match_age_diff_years = ifelse(is.finite(age_diff[idx]), age_diff[idx], NA_real_),
      geob_match_depth_diff_cm = ifelse(is.finite(depth_diff[idx]), depth_diff[idx], NA_real_)
    )
}

nearest_depth_match <- function(layer_df, sparse_df, value_cols, tol_cm = sparse_depth_tolerance_cm) {
  out <- layer_df %>% select(core_clean, depth_in_core_cm, y_bp)
  if (nrow(sparse_df) == 0) {
    for (nm in value_cols) out[[nm]] <- NA_real_
    out$depth_diff_cm <- NA_real_
    out$matched_within_tolerance <- FALSE
    return(out)
  }

nearest_idx <- vapply(
    layer_df$depth_in_core_cm,
    function(d) {
      which.min(abs(sparse_df$depth_in_core_cm - d))
    },
    integer(1)
  )
  depth_diff <- abs(layer_df$depth_in_core_cm - sparse_df$depth_in_core_cm[nearest_idx])
  keep <- is.finite(depth_diff) & depth_diff <= tol_cm

  for (nm in value_cols) {
    vals <- rep(NA_real_, nrow(layer_df))
    vals[keep] <- sparse_df[[nm]][nearest_idx[keep]]
    out[[nm]] <- vals
  }
  out$depth_diff_cm <- depth_diff
  out$matched_within_tolerance <- keep
  out
}

proxy_family <- function(proxy_name) {
  terrig <- c("rb_sr", "al_si", "ti_al", "fe_al", "zr_al", "rb_al")
  prod <- c("ba_ti", "p_ti", "br_ti", "si_ti", "ca_ti")
  redox <- c("mn_fe")
  geob_extra <- c(
    "d13c_per_mil", "toc_percent_weight", "sum_expandable_clays",
    "sum_illites_micas", "sum_montmorillonites_smectites",
    "sum_greenland", "sum_iceland", "ratio_greenland_to_iceland"
  )

  if (proxy_name %in% geob_extra) return("geob25202_extra")
  if (proxy_name %in% redox) return("redox")
  if (proxy_name %in% terrig) return("terrigenous_mineralogical")
  if (proxy_name %in% prod) return("productivity_carbonate")
  if (str_detect(proxy_name, "^(foram_|sst_)")) return("sparse_foram_sst")
  "other"
}

summarise_associations <- function(df, response_col, proxy_cols, analysis_set, response_aggregation) {
  trend_col <- choose_trend(df)

  map_dfr(proxy_cols, function(px) {
    x <- df[[px]]
    y <- df[[response_col]]

    raw <- safe_cor(x, y, method = "spearman")
    det <- list(n = NA_integer_, estimate = NA_real_, p = NA_real_)

    if (!is.na(trend_col)) {
      xr <- safe_ns_residuals(x, df[[trend_col]])
      yr <- safe_ns_residuals(y, df[[trend_col]])
      det <- safe_cor(xr, yr, method = "spearman")
    }

    tibble(
      analysis_set = analysis_set,
      response = response_col,
      response_aggregation = response_aggregation,
      proxy_name = px,
      proxy_family = proxy_family(px),
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

fit_depth_adjusted_models <- function(df, response_cols, proxies, analysis_set = "primary_qc") {
  map_dfr(response_cols, function(resp) {
    map_dfr(proxies, function(px) {
      dat <- df %>%
        select(core_clean, depth_in_core_cm, all_of(resp), all_of(px)) %>%
        filter(is.finite(.data[[resp]]), is.finite(.data[[px]]), is.finite(depth_in_core_cm))

      if (nrow(dat) < min_n_model || dplyr::n_distinct(dat[[px]]) < 4 || dplyr::n_distinct(dat$depth_in_core_cm) < 6) {
        return(tibble(
          analysis_set = analysis_set,
          response = resp,
          proxy_name = px,
          n_model = nrow(dat),
          term = NA_character_,
          estimate = NA_real_,
          std_error = NA_real_,
          p_value = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          adj_r_squared = NA_real_,
          model_note = "not_fit_low_n_or_variation"
        ))
      }

      df_use <- max(1, min(4, dplyr::n_distinct(dat$depth_in_core_cm) - 1))
      dat$proxy_z <- as.numeric(scale(dat[[px]]))
      fit <- tryCatch(
        lm(as.formula(paste0(resp, " ~ splines::ns(depth_in_core_cm, df = ", df_use, ") + proxy_z + core_clean")), data = dat),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(tibble(
          analysis_set = analysis_set,
          response = resp,
          proxy_name = px,
          n_model = nrow(dat),
          term = "proxy_z",
          estimate = NA_real_,
          std_error = NA_real_,
          p_value = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          adj_r_squared = NA_real_,
          model_note = "fit_error"
        ))
      }

      coefs <- summary(fit)$coefficients
      if (!("proxy_z" %in% rownames(coefs))) {
        return(tibble(
          analysis_set = analysis_set,
          response = resp,
          proxy_name = px,
          n_model = nrow(dat),
          term = "proxy_z",
          estimate = NA_real_,
          std_error = NA_real_,
          p_value = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          adj_r_squared = summary(fit)$adj.r.squared,
          model_note = "proxy_term_missing"
        ))
      }

      ci <- tryCatch(confint(fit), error = function(e) NULL)
      tibble(
        analysis_set = analysis_set,
        response = resp,
        proxy_name = px,
        n_model = nrow(dat),
        term = "proxy_z",
        estimate = coefs["proxy_z", "Estimate"],
        std_error = coefs["proxy_z", "Std. Error"],
        p_value = coefs["proxy_z", "Pr(>|t|)"],
        conf_low = if (!is.null(ci)) ci["proxy_z", 1] else NA_real_,
        conf_high = if (!is.null(ci)) ci["proxy_z", 2] else NA_real_,
        adj_r_squared = summary(fit)$adj.r.squared,
        model_note = paste0("depth_ns_df=", df_use, "; core_fixed_effect")
      )
    })
  })
}

# ------------------------------------------------------------------
# Load core datasets
# ------------------------------------------------------------------
damage <- read_tsv(damage_path, show_col_types = FALSE)
classification <- read_tsv(classification_path, show_col_types = FALSE)
metadata <- read_tsv(metadata_path, show_col_types = FALSE)

xrf <- read_tsv(xrf_path, show_col_types = FALSE)
foram <- read_tsv(foram_path, show_col_types = FALSE)
sst <- read_tsv(sst_path, show_col_types = FALSE)
geob <- read_tsv(geob_path, show_col_types = FALSE)

# ------------------------------------------------------------------
# 1) Join diagnostics
# ------------------------------------------------------------------
damage_dup <- damage %>% count(label, name = "n_damage") %>% filter(n_damage > 1)
metadata_dup <- metadata %>% count(label, name = "n_metadata") %>% filter(n_metadata > 1)

join_tbl <- damage %>%
  left_join(
    metadata %>% select(label, core, depth_in_core_cm, y_bp, library_concentration),
    by = "label"
  )

matched_labels <- join_tbl %>% filter(!is.na(core)) %>% distinct(label)
unmatched_labels <- join_tbl %>% filter(is.na(core)) %>% distinct(label)

diag_summary <- tibble(
  diagnostic_type = c(
    "n_damage_labels", "n_metadata_labels", "n_matched_labels", "n_unmatched_labels",
    "n_duplicate_damage_labels", "n_duplicate_metadata_labels"
  ),
  value = c(
    n_distinct(damage$label),
    n_distinct(metadata$label),
    nrow(matched_labels),
    nrow(unmatched_labels),
    nrow(damage_dup),
    nrow(metadata_dup)
  ),
  label = NA_character_,
  n = NA_real_
)

diag_details <- bind_rows(
  unmatched_labels %>% mutate(diagnostic_type = "unmatched_damage_label", value = NA_real_, n = 1),
  damage_dup %>% transmute(diagnostic_type = "duplicate_damage_label", value = NA_real_, label, n = n_damage),
  metadata_dup %>% transmute(diagnostic_type = "duplicate_metadata_label", value = NA_real_, label, n = n_metadata)
) %>%
  select(diagnostic_type, value, label, n)

join_diagnostics <- bind_rows(diag_summary, diag_details)
write_tsv(join_diagnostics, file.path(out_dir, "Ab_join_diagnostics.tsv"))

# ------------------------------------------------------------------
# 2) Build layer-level table (primary + sensitivity)
# ------------------------------------------------------------------
damage_meta <- join_tbl %>%
  mutate(core_clean = normalize_core(core)) %>%
  filter(!is.na(core_clean), is.finite(depth_in_core_cm), is.finite(y_bp))

build_layer <- function(df, analysis_set) {
  df %>%
    group_by(core_clean, depth_in_core_cm, y_bp) %>%
    summarise(
      sample_A_b_median_layer_median = median(sample_A_b_median, na.rm = TRUE),
      sample_A_b_median_layer_mean = mean(sample_A_b_median, na.rm = TRUE),
      sample_A_b_iqr_layer_median = median(sample_A_b_iqr, na.rm = TRUE),
      sample_A_b_iqr_layer_mean = mean(sample_A_b_iqr, na.rm = TRUE),
      n_damage_libraries_per_layer = n(),
      n_rows_sample_median = median(n_rows_sample, na.rm = TRUE),
      n_rows_sample_mean = mean(n_rows_sample, na.rm = TRUE),
      n_rows_sample_min = min(n_rows_sample, na.rm = TRUE),
      n_rows_sample_max = max(n_rows_sample, na.rm = TRUE),
      library_concentration_median = median(library_concentration, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      analysis_set = analysis_set,
      sample_A_b_median_layer_median = ifelse(is.nan(sample_A_b_median_layer_median), NA_real_, sample_A_b_median_layer_median),
      sample_A_b_median_layer_mean = ifelse(is.nan(sample_A_b_median_layer_mean), NA_real_, sample_A_b_median_layer_mean),
      sample_A_b_iqr_layer_median = ifelse(is.nan(sample_A_b_iqr_layer_median), NA_real_, sample_A_b_iqr_layer_median),
      sample_A_b_iqr_layer_mean = ifelse(is.nan(sample_A_b_iqr_layer_mean), NA_real_, sample_A_b_iqr_layer_mean),
      library_concentration_median = ifelse(is.nan(library_concentration_median), NA_real_, library_concentration_median)
    )
}

layer_primary <- damage_meta %>%
  filter(n_rows_sample >= min_rows_primary) %>%
  build_layer("primary_qc")

layer_sensitivity <- damage_meta %>%
  build_layer("sensitivity_all")

layer_summary <- bind_rows(layer_primary, layer_sensitivity) %>%
  arrange(analysis_set, core_clean, depth_in_core_cm)

write_tsv(layer_summary, file.path(out_dir, "Ab_layer_level_summary.tsv"))

# ------------------------------------------------------------------
# 3) Conservative proxy matching
# ------------------------------------------------------------------
xrf <- xrf %>%
  mutate(core_clean = normalize_core(core))

for (i in seq_len(nrow(ratio_definitions))) {
  xrf <- derive_ratio(
    xrf,
    num_col = ratio_definitions$numerator[i],
    den_col = ratio_definitions$denominator[i],
    out_col = ratio_definitions$ratio_name[i]
  )
}

priority_xrf <- ratio_definitions$ratio_name
priority_xrf <- intersect(priority_xrf, names(xrf))

layer_primary_base <- layer_primary
cores_dense <- intersect(unique(layer_primary_base$core_clean), unique(xrf$core_clean))

xrf_layer <- map_dfr(cores_dense, function(cc) {
  layer_core <- layer_primary_base %>% filter(core_clean == cc)
  xrf_core <- xrf %>% filter(core_clean == cc, is.finite(depth_in_core_cm))
  if (nrow(layer_core) == 0) return(tibble())

  agg <- aggregate_dense_to_layers(layer_core, xrf_core, priority_xrf)
  left_join(
    layer_core %>% select(core_clean, depth_in_core_cm, y_bp),
    agg,
    by = "depth_in_core_cm"
  )
})

foram <- foram %>% mutate(core_clean = normalize_core(core))
sst <- sst %>% mutate(core_clean = normalize_core(core))

foram_candidates <- intersect(
  c("polar_planktonic_spp_pct", "sub_polar_planktonic_species_pct", "transitional_planktonic_species_pct", "g_bulloides_pct", "n_pachyderma_pct"),
  names(foram)
)
sst_candidates <- intersect(
  c("sst_uk37_alkenone", "sst_mgca_jonkers_2013", "sst_mgca_kozdon_2009"),
  names(sst)
)

foram_matched <- map_dfr(unique(layer_primary_base$core_clean), function(cc) {
  lc <- layer_primary_base %>% filter(core_clean == cc)
  sc <- foram %>% filter(core_clean == cc, is.finite(depth_in_core_cm))
  if (nrow(lc) == 0) return(tibble())
  nearest_depth_match(lc, sc, foram_candidates, tol_cm = sparse_depth_tolerance_cm) %>%
    rename_with(~paste0("foram_", .x), all_of(foram_candidates)) %>%
    rename(foram_depth_diff_cm = depth_diff_cm, foram_matched_within_tolerance = matched_within_tolerance)
})

sst_matched <- map_dfr(unique(layer_primary_base$core_clean), function(cc) {
  lc <- layer_primary_base %>% filter(core_clean == cc)
  sc <- sst %>% filter(core_clean == cc, is.finite(depth_in_core_cm))
  if (nrow(lc) == 0) return(tibble())
  nearest_depth_match(lc, sc, sst_candidates, tol_cm = sparse_depth_tolerance_cm) %>%
    rename_with(~paste0("sst_", .x), all_of(sst_candidates)) %>%
    rename(sst_depth_diff_cm = depth_diff_cm, sst_matched_within_tolerance = matched_within_tolerance)
})

sparse_quality <- bind_rows(
  foram_matched %>% summarise(
    source = "foram",
    match_fraction = mean(foram_matched_within_tolerance, na.rm = TRUE),
    median_depth_diff_cm = median(foram_depth_diff_cm[foram_matched_within_tolerance], na.rm = TRUE)
  ),
  sst_matched %>% summarise(
    source = "sst",
    match_fraction = mean(sst_matched_within_tolerance, na.rm = TRUE),
    median_depth_diff_cm = median(sst_depth_diff_cm[sst_matched_within_tolerance], na.rm = TRUE)
  )
) %>%
  mutate(
    include_for_analysis = ifelse(
      is.finite(match_fraction) & is.finite(median_depth_diff_cm) &
        match_fraction >= sparse_min_match_fraction &
        median_depth_diff_cm <= sparse_max_median_depth_diff_cm,
      TRUE,
      FALSE
    )
  )

include_foram <- sparse_quality %>% filter(source == "foram") %>% pull(include_for_analysis)
include_sst <- sparse_quality %>% filter(source == "sst") %>% pull(include_for_analysis)
include_foram <- ifelse(length(include_foram) == 0, FALSE, include_foram[1])
include_sst <- ifelse(length(include_sst) == 0, FALSE, include_sst[1])

geob <- geob %>% mutate(core_clean = normalize_core(core))
geob_proxy_cols <- intersect(
  c(
    "d13c_per_mil", "toc_percent_weight", "sum_expandable_clays", "sum_illites_micas",
    "sum_montmorillonites_smectites", "sum_greenland", "sum_iceland", "ratio_greenland_to_iceland"
  ),
  names(geob)
)

geob_layer_match <- layer_primary_base %>%
  filter(core_clean == "GeoB25202") %>%
  rowwise() %>%
  do({
    picked <- nearest_one_row(., geob %>% filter(core_clean == "GeoB25202"))
    if (is.null(picked)) {
      tibble(
        core_clean = .$core_clean,
        depth_in_core_cm = .$depth_in_core_cm,
        y_bp = .$y_bp,
        geob_match_age_diff_years = NA_real_,
        geob_match_depth_diff_cm = NA_real_
      )
    } else {
      tibble(
        core_clean = .$core_clean,
        depth_in_core_cm = .$depth_in_core_cm,
        y_bp = .$y_bp,
        geob_match_age_diff_years = picked$geob_match_age_diff_years,
        geob_match_depth_diff_cm = picked$geob_match_depth_diff_cm
      ) %>%
        bind_cols(picked %>% select(all_of(geob_proxy_cols)))
    }
  }) %>%
  ungroup()

# Build primary matched analysis table
analysis_primary <- layer_primary_base %>%
  left_join(xrf_layer, by = c("core_clean", "depth_in_core_cm", "y_bp")) %>%
  left_join(geob_layer_match, by = c("core_clean", "depth_in_core_cm", "y_bp"))

if (include_foram && nrow(foram_matched) > 0) {
  analysis_primary <- analysis_primary %>%
    left_join(foram_matched %>% select(-foram_matched_within_tolerance), by = c("core_clean", "depth_in_core_cm", "y_bp"))
}
if (include_sst && nrow(sst_matched) > 0) {
  analysis_primary <- analysis_primary %>%
    left_join(sst_matched %>% select(-sst_matched_within_tolerance), by = c("core_clean", "depth_in_core_cm", "y_bp"))
}

# Sensitivity (no QC filter) uses nearest values from primary matched by layer key where available
analysis_sensitivity <- layer_sensitivity %>%
  left_join(
    analysis_primary %>%
      select(-analysis_set) %>%
      distinct(core_clean, depth_in_core_cm, y_bp, .keep_all = TRUE),
    by = c("core_clean", "depth_in_core_cm", "y_bp", "sample_A_b_median_layer_median", "sample_A_b_median_layer_mean", "sample_A_b_iqr_layer_median", "sample_A_b_iqr_layer_mean", "n_damage_libraries_per_layer", "n_rows_sample_median", "n_rows_sample_mean", "n_rows_sample_min", "n_rows_sample_max", "library_concentration_median")
  )

# ------------------------------------------------------------------
# 4) Association summaries
# ------------------------------------------------------------------
priority_proxy_names <- c(
  ratio_definitions$ratio_name,
  geob_proxy_cols,
  if (include_foram) grep("^foram_", names(analysis_primary), value = TRUE) else character(0),
  if (include_sst) grep("^sst_", names(analysis_primary), value = TRUE) else character(0)
)
priority_proxy_names <- intersect(unique(priority_proxy_names), names(analysis_primary))

assoc_primary <- bind_rows(
  summarise_associations(
    analysis_primary,
    response_col = "sample_A_b_median_layer_median",
    proxy_cols = priority_proxy_names,
    analysis_set = "primary_qc",
    response_aggregation = "median"
  ),
  summarise_associations(
    analysis_primary,
    response_col = "sample_A_b_iqr_layer_median",
    proxy_cols = priority_proxy_names,
    analysis_set = "primary_qc",
    response_aggregation = "median"
  ),
  summarise_associations(
    analysis_primary,
    response_col = "sample_A_b_median_layer_mean",
    proxy_cols = priority_proxy_names,
    analysis_set = "primary_qc",
    response_aggregation = "mean_sensitivity"
  ),
  summarise_associations(
    analysis_primary,
    response_col = "sample_A_b_iqr_layer_mean",
    proxy_cols = priority_proxy_names,
    analysis_set = "primary_qc",
    response_aggregation = "mean_sensitivity"
  )
)

assoc_sensitivity <- bind_rows(
  summarise_associations(
    analysis_sensitivity,
    response_col = "sample_A_b_median_layer_median",
    proxy_cols = priority_proxy_names,
    analysis_set = "sensitivity_all",
    response_aggregation = "median"
  ),
  summarise_associations(
    analysis_sensitivity,
    response_col = "sample_A_b_iqr_layer_median",
    proxy_cols = priority_proxy_names,
    analysis_set = "sensitivity_all",
    response_aggregation = "median"
  )
)

assoc_summary <- bind_rows(assoc_primary, assoc_sensitivity) %>%
  mutate(
    q_raw_bh = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  arrange(response, analysis_set, desc(abs(spearman_detrended)), p_detrended)

write_tsv(assoc_summary, file.path(out_dir, "Ab_proxy_association_summary.tsv"))

# ------------------------------------------------------------------
# 5) Selected depth-adjusted spline models (high priority)
# ------------------------------------------------------------------
high_priority_model_proxies <- intersect(
  c("rb_sr", "al_si", "ti_al", "fe_al", "ba_ti", "p_ti", "br_ti", "si_ti", "ca_ti", "mn_fe", "ratio_greenland_to_iceland", "toc_percent_weight", "sum_expandable_clays"),
  names(analysis_primary)
)

model_summary <- fit_depth_adjusted_models(
  analysis_primary,
  response_cols = c("sample_A_b_median_layer_median", "sample_A_b_iqr_layer_median"),
  proxies = high_priority_model_proxies,
  analysis_set = "primary_qc"
)
write_tsv(model_summary, file.path(out_dir, "Ab_selected_model_summary.tsv"))

# ------------------------------------------------------------------
# 6) GeoB25202 extra-proxy focused summary
# ------------------------------------------------------------------
geob_only <- analysis_primary %>% filter(core_clean == "GeoB25202")
geob_extra_summary <- bind_rows(
  summarise_associations(
    geob_only,
    response_col = "sample_A_b_median_layer_median",
    proxy_cols = intersect(geob_proxy_cols, names(geob_only)),
    analysis_set = "geob25202_only",
    response_aggregation = "median"
  ),
  summarise_associations(
    geob_only,
    response_col = "sample_A_b_iqr_layer_median",
    proxy_cols = intersect(geob_proxy_cols, names(geob_only)),
    analysis_set = "geob25202_only",
    response_aggregation = "median"
  )
) %>%
  mutate(
    q_raw_bh = ifelse(is.na(p_raw), NA_real_, p.adjust(p_raw, method = "BH")),
    q_detrended_bh = ifelse(is.na(p_detrended), NA_real_, p.adjust(p_detrended, method = "BH"))
  ) %>%
  arrange(response, desc(abs(spearman_detrended)), p_detrended)

write_tsv(geob_extra_summary, file.path(out_dir, "Ab_geob25202_extra_proxy_summary.tsv"))

# ------------------------------------------------------------------
# 7) Plots
# ------------------------------------------------------------------
context_df <- analysis_primary %>%
  select(core_clean, depth_in_core_cm, y_bp, sample_A_b_median_layer_median, sample_A_b_iqr_layer_median) %>%
  pivot_longer(
    cols = c(sample_A_b_median_layer_median, sample_A_b_iqr_layer_median),
    names_to = "response",
    values_to = "value"
  )

p_context <- ggplot(context_df, aes(x = y_bp, y = value, color = core_clean)) +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_line(alpha = 0.7) +
  scale_x_reverse() +
  facet_wrap(~response, scales = "free_y", ncol = 1) +
  labs(
    x = "Age (y BP)",
    y = "Layer-level A_b metric",
    color = "Core",
    title = "A_b context across cores (primary QC layers)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "Ab_multicore_context_over_age.png"),
  plot = p_context,
  width = 10,
  height = 7,
  dpi = 300
)

corr_plot_df <- assoc_summary %>%
  filter(
    analysis_set == "primary_qc",
    response == "sample_A_b_median_layer_median",
    response_aggregation == "median",
    is.finite(spearman_raw)
  ) %>%
  select(proxy_name, proxy_family, n_raw, spearman_raw, spearman_detrended) %>%
  pivot_longer(cols = c(spearman_raw, spearman_detrended), names_to = "metric", values_to = "rho")

p_corr <- ggplot(corr_plot_df, aes(x = rho, y = reorder(proxy_name, rho), color = proxy_family, shape = metric)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(size = 2.3, alpha = 0.9) +
  labs(
    x = "Spearman correlation",
    y = "Proxy",
    color = "Proxy family",
    shape = "Association type",
    title = "A_b median proxy associations (raw vs detrended)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "Ab_median_proxy_correlation_summary.png"),
  plot = p_corr,
  width = 10,
  height = 8,
  dpi = 300
)

if (nrow(geob_extra_summary) > 0) {
  p_geob <- geob_extra_summary %>%
    filter(response == "sample_A_b_median_layer_median") %>%
    pivot_longer(cols = c(spearman_raw, spearman_detrended), names_to = "metric", values_to = "rho") %>%
    ggplot(aes(x = rho, y = reorder(proxy_name, rho), shape = metric)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(size = 2.5, color = "#0D47A1") +
    labs(
      x = "Spearman correlation",
      y = "GeoB25202 extra proxy",
      shape = "Association type",
      title = "GeoB25202 extra proxy summary for A_b median"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(plot_dir, "Ab_geob25202_extra_proxy_summary.png"),
    plot = p_geob,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ------------------------------------------------------------------
# 8) Notes (guardrails + assumptions)
# ------------------------------------------------------------------
top_primary_median <- assoc_summary %>%
  filter(
    analysis_set == "primary_qc",
    response == "sample_A_b_median_layer_median",
    response_aggregation == "median",
    is.finite(spearman_detrended)
  ) %>%
  arrange(desc(abs(spearman_detrended))) %>%
  slice_head(n = 8)

notes <- c(
  "A_b interpretation guardrails:",
  "- sample_A_b_* are damage-related metrics, not direct ecological abundance proxies.",
  "- observed associations with paleo/sediment proxies may reflect preservation environment, sediment matrix, age/depth structure, or authenticity-related processes.",
  "- these results are observational and non-causal.",
  "- higher A_b should not be conflated with improved ecological interpretability without additional authentication context.",
  "- XRF input for this run is data/combined_xrf_not_normalized.tsv and is interpreted only through element ratios.",
  "- Raw non-normalized XRF element values were excluded from priority/model summaries.",
  "",
  "Implementation assumptions and conservative choices:",
  paste0("- Primary QC filter used: n_rows_sample >= ", min_rows_primary, "."),
  "- Layer collapse key: core_clean + depth_in_core_cm + y_bp.",
  "- Technical replicate collapse for A_b used median as primary and mean as sensitivity.",
  paste0("- Core names normalized by stripping replicate suffix _R1/_R2.")
)

if (!file.exists(classification_path_primary) && file.exists(classification_path_fallback)) {
  notes <- c(
    notes,
    "- sample-classification-stats.tsv was loaded from exploration/ fallback path."
  )
}

if (nrow(sparse_quality) > 0) {
  notes <- c(
    notes,
    "",
    "Sparse match quality diagnostics:",
    apply(sparse_quality, 1, function(z) {
      paste0(
        "- ", z[["source"]], ": match_fraction=", sprintf("%.3f", as.numeric(z[["match_fraction"]])),
        ", median_depth_diff_cm=", sprintf("%.3f", as.numeric(z[["median_depth_diff_cm"]])),
        ", included=", z[["include_for_analysis"]]
      )
    })
  )
}

if (nrow(top_primary_median) > 0) {
  notes <- c(
    notes,
    "",
    "Top detrended associations for A_b median (primary QC; ranked by |rho|):",
    apply(top_primary_median, 1, function(z) {
      paste0(
        "- ", z[["proxy_name"]], " (", z[["proxy_family"]], "): rho_det=",
        sprintf("%.3f", as.numeric(z[["spearman_detrended"]])),
        ", p_det=", sprintf("%.3g", as.numeric(z[["p_detrended"]])),
        ", n=", z[["n_detrended"]]
      )
    })
  )
}

writeLines(notes, con = file.path(out_dir, "Ab_results_notes.txt"))

message("Done. Outputs written to: ", out_dir)
message("Plots written to: ", plot_dir)
