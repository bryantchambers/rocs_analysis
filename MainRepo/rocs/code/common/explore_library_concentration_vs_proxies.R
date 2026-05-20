#!/usr/bin/env Rscript

# Exploratory analysis only:
# - tests associations between library_concentration and paleo proxies
# - does not imply causation
# - pooled results may be confounded by age structure and core differences

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

# -------------------------------------------------------------------
# Parameters (kept explicit for reproducibility)
# -------------------------------------------------------------------
age_tolerance_years <- 500
min_matched_samples <- 10
min_matched_samples_within_core <- 6

out_dir <- file.path("analysis", "library_concentration_proxy_exploration")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path("data", "metadata_v5.tsv")
xrf_path <- file.path("data", "combined_xrf_geochemistry_curated.csv")
foram_path <- file.path("data", "combined_foraminifera_geochem.tsv")
sst_path <- file.path("data", "combined_sst_proxies_separate_columns.tsv")
geob_path <- file.path("data", "geob25202_clean_proxies.tsv")

required_files <- c(metadata_path, xrf_path, foram_path, sst_path, geob_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input files (run from repo root): ",
    paste(missing_files, collapse = ", ")
  )
}

normalize_core <- function(x) {
  str_replace(x, "_R[12]$", "")
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x_ok <- x[ok]
  y_ok <- y[ok]
  if (length(x_ok) < 3 || dplyr::n_distinct(x_ok) < 2 || dplyr::n_distinct(y_ok) < 2) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  out <- suppressWarnings(cor.test(x_ok, y_ok, method = method, exact = FALSE))
  list(estimate = unname(out$estimate), p.value = out$p.value)
}

make_proxy_long <- function(df, source_name, y_bp_col = "y_bp", exclude_numeric = character()) {
  if (!("core" %in% names(df))) stop(source_name, ": missing 'core' column")
  if (!(y_bp_col %in% names(df))) stop(source_name, ": missing age column '", y_bp_col, "'")

  df <- df %>%
    mutate(
      core_clean = normalize_core(.data$core),
      y_bp_proxy = as.numeric(.data[[y_bp_col]])
    )

  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  proxy_cols <- setdiff(numeric_cols, c(y_bp_col, "y_bp_proxy", exclude_numeric))

  df %>%
    select(core_clean, y_bp_proxy, all_of(proxy_cols)) %>%
    pivot_longer(
      cols = all_of(proxy_cols),
      names_to = "proxy_name",
      values_to = "proxy_value"
    ) %>%
    mutate(source_dataset = source_name)
}

match_proxy_group <- function(sample_df, proxy_group_df, tolerance_years) {
  core_value <- unique(proxy_group_df$core_clean)
  if (length(core_value) != 1) {
    return(tibble())
  }

  sample_core <- sample_df %>% filter(core_clean == core_value)
  proxy_group_df <- proxy_group_df %>%
    filter(is.finite(y_bp_proxy), is.finite(proxy_value)) %>%
    arrange(y_bp_proxy)

  if (nrow(sample_core) == 0 || nrow(proxy_group_df) == 0) {
    return(tibble())
  }

  proxy_ages <- proxy_group_df$y_bp_proxy
  nearest_idx <- vapply(
    sample_core$y_bp,
    function(age_sample) {
      if (!is.finite(age_sample) || length(proxy_ages) == 0) return(NA_integer_)
      which.min(abs(proxy_ages - age_sample))
    },
    integer(1)
  )

  nearest_age <- proxy_ages[nearest_idx]
  age_diff <- abs(sample_core$y_bp - nearest_age)
  keep <- is.finite(age_diff) & age_diff <= tolerance_years

  if (!any(keep)) {
    return(tibble())
  }

  proxy_pick <- proxy_group_df[nearest_idx[keep], c("source_dataset", "proxy_name", "y_bp_proxy", "proxy_value")]

  bind_cols(
    sample_core[keep, ],
    proxy_pick,
    tibble(age_diff = age_diff[keep])
  )
}

summarise_correlations <- function(df, min_n = 10) {
  out <- df %>%
    group_by(source_dataset, proxy_name) %>%
    reframe({
      x <- proxy_value
      y <- library_concentration_log10
      y_raw <- library_concentration_mean

      n_use <- sum(is.finite(x) & is.finite(y))
      if (n_use < min_n) {
        tibble(
          n_matched = n_use,
          spearman_rho = NA_real_,
          spearman_p = NA_real_,
          kendall_tau = NA_real_,
          kendall_p = NA_real_,
          median_age_diff = median(age_diff, na.rm = TRUE),
          max_age_diff = max(age_diff, na.rm = TRUE),
          library_vs_age_spearman = safe_cor_test(y_raw, y_bp, method = "spearman")$estimate,
          proxy_vs_age_spearman = safe_cor_test(x, y_bp, method = "spearman")$estimate
        )
      } else {
        sp <- safe_cor_test(x, y, method = "spearman")
        kd <- safe_cor_test(x, y, method = "kendall")
        tibble(
          n_matched = n_use,
          spearman_rho = sp$estimate,
          spearman_p = sp$p.value,
          kendall_tau = kd$estimate,
          kendall_p = kd$p.value,
          median_age_diff = median(age_diff, na.rm = TRUE),
          max_age_diff = max(age_diff, na.rm = TRUE),
          library_vs_age_spearman = safe_cor_test(y_raw, y_bp, method = "spearman")$estimate,
          proxy_vs_age_spearman = safe_cor_test(x, y_bp, method = "spearman")$estimate
        )
      }
    }) %>%
    ungroup() %>%
    mutate(
      q_value_bh = if_else(
        is.na(spearman_p),
        NA_real_,
        p.adjust(spearman_p, method = "BH")
      )
    ) %>%
    arrange(desc(abs(spearman_rho)), spearman_p)

  out
}

# -------------------------------------------------------------------
# 1) Metadata: collapse technical replicates to one sample-layer row
# -------------------------------------------------------------------
metadata_raw <- read_tsv(metadata_path, show_col_types = FALSE)

metadata_collapsed <- metadata_raw %>%
  mutate(core_clean = normalize_core(core)) %>%
  filter(!is.na(core_clean), !is.na(y_bp), !is.na(depth_in_core_cm)) %>%
  group_by(core_clean, y_bp, depth_in_core_cm) %>%
  summarise(
    n_replicates = n(),
    n_replicates_library_nonmissing = sum(!is.na(library_concentration)),
    library_concentration_mean = mean(library_concentration, na.rm = TRUE),
    library_concentration_median = median(library_concentration, na.rm = TRUE),
    temp_mean = mean(temp, na.rm = TRUE),
    mis_mean = mean(mis, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    library_concentration_mean = if_else(is.nan(library_concentration_mean), NA_real_, library_concentration_mean),
    library_concentration_median = if_else(is.nan(library_concentration_median), NA_real_, library_concentration_median),
    temp_mean = if_else(is.nan(temp_mean), NA_real_, temp_mean),
    mis_mean = if_else(is.nan(mis_mean), NA_real_, mis_mean)
  )

min_positive_lib <- suppressWarnings(min(metadata_collapsed$library_concentration_mean[metadata_collapsed$library_concentration_mean > 0], na.rm = TRUE))
if (!is.finite(min_positive_lib)) min_positive_lib <- 1e-6
small_value <- min_positive_lib / 2

metadata_collapsed <- metadata_collapsed %>%
  mutate(
    library_concentration_log10 = log10(library_concentration_mean + small_value),
    sample_layer_id = paste(core_clean, y_bp, depth_in_core_cm, sep = "|")
  )

# -------------------------------------------------------------------
# 2) Build long-format proxy table from all proxy sources
# -------------------------------------------------------------------
xrf <- read_csv(xrf_path, show_col_types = FALSE)
foram <- read_tsv(foram_path, show_col_types = FALSE)
geob <- read_tsv(geob_path, show_col_types = FALSE)
sst_raw <- read_tsv(sst_path, show_col_types = FALSE)

# SST file may not contain y_bp; map by core + depth when needed
if (!("y_bp" %in% names(sst_raw))) {
  sst <- sst_raw %>%
    left_join(
      geob %>%
        select(core, depth_in_core_cm, y_bp) %>%
        distinct(),
      by = c("core", "depth_in_core_cm")
    )
} else {
  sst <- sst_raw
}

proxy_long <- bind_rows(
  make_proxy_long(
    xrf,
    source_name = "combined_xrf_geochemistry_curated.csv",
    y_bp_col = "y_bp",
    exclude_numeric = c("depth_in_core_cm")
  ),
  make_proxy_long(
    foram,
    source_name = "combined_foraminifera_geochem.tsv",
    y_bp_col = "y_bp",
    exclude_numeric = c("depth_in_core_cm")
  ),
  make_proxy_long(
    sst,
    source_name = "combined_sst_proxies_separate_columns.tsv",
    y_bp_col = "y_bp",
    exclude_numeric = c("depth_in_core_cm")
  ),
  make_proxy_long(
    geob,
    source_name = "geob25202_clean_proxies.tsv",
    y_bp_col = "y_bp",
    exclude_numeric = c("depth_in_core_cm")
  ),
  metadata_collapsed %>%
    transmute(
      core_clean,
      y_bp_proxy = y_bp,
      source_dataset = "metadata_v5.tsv",
      temp = temp_mean,
      mis = mis_mean
    ) %>%
    pivot_longer(cols = c(temp, mis), names_to = "proxy_name", values_to = "proxy_value")
) %>%
  filter(is.finite(y_bp_proxy), is.finite(proxy_value), !is.na(core_clean))

# -------------------------------------------------------------------
# 3) Nearest-age matching within core and conservative tolerance
# -------------------------------------------------------------------
sample_layers <- metadata_collapsed %>%
  filter(is.finite(library_concentration_mean), is.finite(y_bp)) %>%
  select(sample_layer_id, core_clean, y_bp, depth_in_core_cm, n_replicates,
         library_concentration_mean, library_concentration_log10)

proxy_groups <- proxy_long %>%
  group_by(source_dataset, proxy_name, core_clean) %>%
  group_split(.keep = TRUE)

matched <- map_dfr(
  proxy_groups,
  ~match_proxy_group(
    sample_df = sample_layers,
    proxy_group_df = .x,
    tolerance_years = age_tolerance_years
  )
)

if (nrow(matched) == 0) {
  stop("No sample-proxy matches found within tolerance = ", age_tolerance_years, " years.")
}

# -------------------------------------------------------------------
# 4) Correlation summaries (overall + within-core)
# -------------------------------------------------------------------
cor_overall <- summarise_correlations(matched, min_n = min_matched_samples)

cor_within_core <- matched %>%
  group_by(core_clean, source_dataset, proxy_name) %>%
  reframe({
    x <- proxy_value
    y <- library_concentration_log10
    n_use <- sum(is.finite(x) & is.finite(y))

    if (n_use < min_matched_samples_within_core) {
      tibble(
        n_matched = n_use,
        spearman_rho = NA_real_,
        spearman_p = NA_real_,
        kendall_tau = NA_real_,
        kendall_p = NA_real_,
        median_age_diff = median(age_diff, na.rm = TRUE),
        max_age_diff = max(age_diff, na.rm = TRUE)
      )
    } else {
      sp <- safe_cor_test(x, y, method = "spearman")
      kd <- safe_cor_test(x, y, method = "kendall")
      tibble(
        n_matched = n_use,
        spearman_rho = sp$estimate,
        spearman_p = sp$p.value,
        kendall_tau = kd$estimate,
        kendall_p = kd$p.value,
        median_age_diff = median(age_diff, na.rm = TRUE),
        max_age_diff = max(age_diff, na.rm = TRUE)
      )
    }
  }) %>%
  ungroup() %>%
  group_by(core_clean) %>%
  mutate(q_value_bh = if_else(is.na(spearman_p), NA_real_, p.adjust(spearman_p, method = "BH"))) %>%
  ungroup() %>%
  arrange(core_clean, desc(abs(spearman_rho)), spearman_p)

# -------------------------------------------------------------------
# 5) Output tables
# -------------------------------------------------------------------
write_tsv(
  metadata_collapsed,
  file.path(out_dir, "metadata_collapsed_library_layers.tsv")
)

write_tsv(
  matched %>%
    select(sample_layer_id, core_clean, y_bp, depth_in_core_cm,
           n_replicates, library_concentration_mean, library_concentration_log10,
           source_dataset, proxy_name, y_bp_proxy, age_diff, proxy_value),
  file.path(out_dir, "matched_sample_proxy_pairs.tsv")
)

write_tsv(cor_overall, file.path(out_dir, "correlations_overall.tsv"))
write_tsv(cor_within_core, file.path(out_dir, "correlations_within_core.tsv"))

# -------------------------------------------------------------------
# 6) Exploratory plots
# -------------------------------------------------------------------
p_hist <- metadata_collapsed %>%
  filter(is.finite(library_concentration_mean)) %>%
  ggplot(aes(x = library_concentration_mean)) +
  geom_histogram(bins = 40, fill = "#2C7FB8", color = "white", alpha = 0.9) +
  geom_density(aes(y = after_stat(count)), color = "#253494", linewidth = 0.7, na.rm = TRUE) +
  labs(
    title = "Distribution of library concentration (collapsed technical replicates)",
    subtitle = "Exploratory context only",
    x = "Library concentration (nM)",
    y = "Count"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "library_concentration_distribution.png"),
  plot = p_hist,
  width = 8,
  height = 5,
  dpi = 300
)

top_tested <- cor_overall %>%
  filter(is.finite(spearman_rho), n_matched >= min_matched_samples) %>%
  arrange(desc(spearman_rho))

top_pos <- head(top_tested, 3)
top_neg <- head(arrange(top_tested, spearman_rho), 3)
top_scatter <- bind_rows(top_pos, top_neg) %>%
  distinct(source_dataset, proxy_name, .keep_all = TRUE) %>%
  mutate(
    proxy_label = paste0(
      proxy_name,
      "\n", source_dataset,
      "\n", "rho=", formatC(spearman_rho, format = "f", digits = 2),
      ", p=", formatC(spearman_p, format = "e", digits = 2),
      ", n=", n_matched
    )
  )

if (nrow(top_scatter) > 0) {
  scatter_df <- matched %>%
    inner_join(top_scatter %>% select(source_dataset, proxy_name, proxy_label),
               by = c("source_dataset", "proxy_name"))

  p_scatter <- scatter_df %>%
    ggplot(aes(x = proxy_value, y = library_concentration_log10, color = core_clean)) +
    geom_point(alpha = 0.75, size = 1.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    facet_wrap(~proxy_label, scales = "free_x") +
    labs(
      title = "Top positive and negative exploratory proxy associations",
      subtitle = "Spearman ranking on pooled matched layers (exploratory only)",
      x = "Proxy value",
      y = "log10(library concentration + offset)",
      color = "Core"
    ) +
    theme_bw() +
    theme(strip.text = element_text(size = 8))

  ggsave(
    filename = file.path(plot_dir, "top_positive_negative_proxy_scatter.png"),
    plot = p_scatter,
    width = 12,
    height = 7,
    dpi = 300
  )
}

confounding_df <- cor_overall %>%
  filter(is.finite(spearman_rho), n_matched >= min_matched_samples)

if (nrow(confounding_df) > 0) {
  label_df <- confounding_df %>%
    arrange(desc(abs(spearman_rho))) %>%
    slice_head(n = min(12, nrow(confounding_df)))

  p_confound <- confounding_df %>%
    ggplot(aes(
      x = spearman_rho,
      y = proxy_vs_age_spearman,
      size = n_matched,
      color = source_dataset
    )) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_point(alpha = 0.75) +
    geom_text(
      data = label_df,
      aes(label = proxy_name),
      size = 2.6,
      vjust = -0.6,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    labs(
      title = "Age-structure diagnostic for exploratory proxy associations",
      subtitle = "X: proxy vs log10(library concentration), Y: proxy vs age (Spearman)",
      x = "Spearman rho (proxy vs log10 library concentration)",
      y = "Spearman rho (proxy vs age)",
      color = "Source",
      size = "Matched n"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(plot_dir, "age_structure_diagnostic.png"),
    plot = p_confound,
    width = 10,
    height = 7,
    dpi = 300
  )
}

message("Done. Exploratory outputs written to: ", out_dir)
message("Conservative caveat: pooled associations may reflect age/core structure and age-matching uncertainty.")
