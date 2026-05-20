library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(tibble)

ROCS_ROOT <- Sys.getenv("ROCS_ROOT", "/projects/caeg/people/ngm902/apps/repos/rocs")
ROCS_ROOT <- normalizePath(ROCS_ROOT, mustWork = FALSE)

params <- list(
  model_version = "damage-baseline-proxy-association-v1",
  depth_match_tolerance_cm = 1.0,
  min_pairs_for_test = 20,
  top_hits_per_response = 20,
  random_seed = 20260325
)

set.seed(params$random_seed)

input_sample_baselines <- file.path(ROCS_ROOT, "results/microbial/damage/damage-classification-depositional/sample-baselines.tsv")
input_metadata <- file.path(ROCS_ROOT, "data/metadata_v5.tsv")
input_forams <- file.path(ROCS_ROOT, "data/combined_foraminifera_geochem.tsv")
input_sst <- file.path(ROCS_ROOT, "data/combined_sst_proxies_separate_columns.tsv")
input_xrf <- file.path(ROCS_ROOT, "data/combined_xrf_geochemistry_curated.csv")
input_geob <- file.path(ROCS_ROOT, "data/geob25202_clean_proxies.tsv")

out_dir <- file.path(ROCS_ROOT, "results/microbial/damage/damage-baseline-proxy-associations")
plots_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")

log_msg <- function(...) {
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = log_file, append = TRUE)
}

read_table_auto <- function(path) {
  if (str_detect(path, "\\.csv$")) {
    read_csv(path, show_col_types = FALSE)
  } else {
    read_tsv(path, show_col_types = FALSE)
  }
}

as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

prefix_proxy_columns <- function(df, id_cols, prefix) {
  proxy_cols <- setdiff(names(df), id_cols)
  names(df)[match(proxy_cols, names(df))] <- paste0(prefix, proxy_cols)
  df
}

nearest_join_by_depth <- function(target_df, proxy_df, proxy_name, tolerance_cm = 1) {
  id_cols <- c("core_base", "depth_in_core_cm")
  if (!all(id_cols %in% names(target_df))) {
    stop("target_df must contain core_base and depth_in_core_cm")
  }
  if (!all(id_cols %in% names(proxy_df))) {
    stop("proxy_df must contain core_base and depth_in_core_cm")
  }

  proxy_vars <- setdiff(names(proxy_df), id_cols)
  out <- target_df
  for (col in proxy_vars) {
    out[[col]] <- NA_real_
  }
  dist_col <- paste0(proxy_name, "_match_dist_cm")
  out[[dist_col]] <- NA_real_

  common_cores <- intersect(unique(target_df$core_base), unique(proxy_df$core_base))
  if (length(common_cores) == 0) {
    return(out)
  }

  for (cc in common_cores) {
    i_target <- which(target_df$core_base == cc)
    i_proxy <- which(proxy_df$core_base == cc)
    proxy_depth <- proxy_df$depth_in_core_cm[i_proxy]

    for (ii in i_target) {
      d <- abs(proxy_depth - target_df$depth_in_core_cm[ii])
      if (all(is.na(d))) next
      j_local <- which.min(d)
      if (length(j_local) == 0 || is.na(d[j_local]) || d[j_local] > tolerance_cm) next

      j_proxy <- i_proxy[j_local]
      out[ii, proxy_vars] <- proxy_df[j_proxy, proxy_vars]
      out[[dist_col]][ii] <- d[j_local]
    }
  }

  out
}

compute_association <- function(df, response, predictor, min_pairs = 20) {
  dat <- df %>%
    select(core_base, y_bp, all_of(response), all_of(predictor)) %>%
    mutate(
      y = as_numeric_safe(.data[[response]]),
      x = as_numeric_safe(.data[[predictor]])
    ) %>%
    filter(!is.na(y), !is.na(x), is.finite(y), is.finite(x))

  n_pairs <- nrow(dat)
  n_cores <- n_distinct(dat$core_base)
  if (n_pairs < min_pairs || length(unique(dat$x)) < 3 || length(unique(dat$y)) < 3) {
    return(tibble(
      n_pairs = n_pairs,
      n_cores = n_cores,
      spearman_rho = NA_real_,
      spearman_p = NA_real_,
      lm_beta = NA_real_,
      lm_p = NA_real_,
      lm_adj_r2 = NA_real_
    ))
  }

  spearman_test <- suppressWarnings(cor.test(dat$x, dat$y, method = "spearman", exact = FALSE))

  if (n_cores >= 2) {
    fit <- lm(scale(y) ~ scale(x) + poly(y_bp, 2, raw = TRUE) + core_base, data = dat)
  } else {
    fit <- lm(scale(y) ~ scale(x) + poly(y_bp, 2, raw = TRUE), data = dat)
  }

  fit_sum <- summary(fit)
  coef_names <- rownames(coef(fit_sum))
  x_row <- which(coef_names == "scale(x)")

  lm_beta <- if (length(x_row) == 1) coef(fit_sum)[x_row, "Estimate"] else NA_real_
  lm_p <- if (length(x_row) == 1) coef(fit_sum)[x_row, "Pr(>|t|)"] else NA_real_

  tibble(
    n_pairs = n_pairs,
    n_cores = n_cores,
    spearman_rho = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value,
    lm_beta = lm_beta,
    lm_p = lm_p,
    lm_adj_r2 = unname(fit_sum$adj.r.squared)
  )
}

log_msg("Reading sample baselines and metadata")
sample_baselines <- read_tsv(input_sample_baselines, show_col_types = FALSE)
metadata <- read_tsv(input_metadata, show_col_types = FALSE)

required_baseline_cols <- c("label", "sample_A_b_median", "sample_A_b_iqr", "sample_local_fixed_rate", "sample_Zfit_median", "sample_rho_median")
missing_baseline_cols <- setdiff(required_baseline_cols, names(sample_baselines))
if (length(missing_baseline_cols) > 0) {
  stop("Missing required columns in sample-baselines: ", paste(missing_baseline_cols, collapse = ", "))
}

required_metadata_cols <- c("label", "core", "depth_in_core_cm", "y_bp")
missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata))
if (length(missing_metadata_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_metadata_cols, collapse = ", "))
}

damage_with_meta <- sample_baselines %>%
  left_join(metadata %>% select(label, core, depth_in_core_cm, y_bp), by = "label") %>%
  mutate(
    core_base = str_remove(core, "_R[12]$"),
    depth_in_core_cm = as_numeric_safe(depth_in_core_cm),
    y_bp = as_numeric_safe(y_bp)
  )

unmatched_labels <- sum(is.na(damage_with_meta$core))
log_msg("Rows in sample baselines: ", nrow(sample_baselines))
log_msg("Rows without metadata match: ", unmatched_labels)

damage_agg <- damage_with_meta %>%
  filter(!is.na(core_base), !is.na(depth_in_core_cm)) %>%
  group_by(core_base, depth_in_core_cm) %>%
  summarise(
    y_bp = median(y_bp, na.rm = TRUE),
    sample_A_b_median = median(sample_A_b_median, na.rm = TRUE),
    sample_A_b_iqr = median(sample_A_b_iqr, na.rm = TRUE),
    sample_local_fixed_rate = median(sample_local_fixed_rate, na.rm = TRUE),
    sample_Zfit_median = median(sample_Zfit_median, na.rm = TRUE),
    sample_rho_median = median(sample_rho_median, na.rm = TRUE),
    n_labels_collapsed = n(),
    .groups = "drop"
  )

log_msg("Collapsed damage rows (core-depth): ", nrow(damage_agg))

log_msg("Reading proxy tables")
forams <- read_table_auto(input_forams) %>%
  mutate(
    core_base = as.character(core),
    depth_in_core_cm = as_numeric_safe(depth_in_core_cm)
  ) %>%
  select(-any_of(c("core", "y_bp")))

sst <- read_table_auto(input_sst) %>%
  mutate(
    core_base = as.character(core),
    depth_in_core_cm = as_numeric_safe(depth_in_core_cm)
  ) %>%
  select(-any_of(c("core")))

xrf <- read_table_auto(input_xrf) %>%
  mutate(
    core_base = as.character(core),
    depth_in_core_cm = as_numeric_safe(depth_in_core_cm)
  ) %>%
  select(-any_of(c("core", "y_bp")))

geob <- read_table_auto(input_geob) %>%
  mutate(
    core_base = as.character(core),
    depth_in_core_cm = as_numeric_safe(depth_in_core_cm)
  ) %>%
  select(-any_of(c("core", "y_bp")))

forams <- prefix_proxy_columns(forams, c("core_base", "depth_in_core_cm"), "foram_")
sst <- prefix_proxy_columns(sst, c("core_base", "depth_in_core_cm"), "sst_")
xrf <- prefix_proxy_columns(xrf, c("core_base", "depth_in_core_cm"), "xrf_")
geob <- prefix_proxy_columns(geob, c("core_base", "depth_in_core_cm"), "geob_")

proxy_joined <- damage_agg %>%
  nearest_join_by_depth(forams, proxy_name = "foram", tolerance_cm = params$depth_match_tolerance_cm) %>%
  nearest_join_by_depth(sst, proxy_name = "sst", tolerance_cm = params$depth_match_tolerance_cm) %>%
  nearest_join_by_depth(xrf, proxy_name = "xrf", tolerance_cm = params$depth_match_tolerance_cm) %>%
  nearest_join_by_depth(geob, proxy_name = "geob", tolerance_cm = params$depth_match_tolerance_cm)

response_vars <- c(
  "sample_A_b_median",
  "sample_A_b_iqr",
  "sample_local_fixed_rate",
  "sample_Zfit_median",
  "sample_rho_median"
)

exclude_cols <- c("core_base", "depth_in_core_cm", "y_bp", response_vars, "n_labels_collapsed")
predictor_vars <- setdiff(names(proxy_joined), exclude_cols)
predictor_vars <- predictor_vars[!str_detect(predictor_vars, "_match_dist_cm$")]

log_msg("Number of candidate predictors: ", length(predictor_vars))

association_results <- crossing(response = response_vars, predictor = predictor_vars) %>%
  mutate(res = map2(response, predictor, ~ compute_association(proxy_joined, .x, .y, params$min_pairs_for_test))) %>%
  unnest(res) %>%
  mutate(
    source = case_when(
      str_starts(predictor, "foram_") ~ "foraminifera_geochem",
      str_starts(predictor, "sst_") ~ "sst",
      str_starts(predictor, "xrf_") ~ "xrf_geochem",
      str_starts(predictor, "geob_") ~ "geob25202_clean",
      TRUE ~ "unknown"
    )
  ) %>%
  group_by(response) %>%
  mutate(
    spearman_fdr = p.adjust(spearman_p, method = "BH"),
    lm_fdr = p.adjust(lm_p, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(response, lm_fdr, spearman_fdr, desc(abs(spearman_rho)))

top_hits <- association_results %>%
  filter(!is.na(spearman_rho)) %>%
  group_by(response) %>%
  arrange(lm_fdr, spearman_fdr, desc(abs(spearman_rho)), .by_group = TRUE) %>%
  slice_head(n = params$top_hits_per_response) %>%
  ungroup()

log_msg("Creating plots for top hits")
walk(unique(top_hits$response), function(resp) {
  top_preds <- top_hits %>%
    filter(response == resp) %>%
    slice_head(n = 6) %>%
    pull(predictor)

  if (length(top_preds) == 0) return(NULL)

  plot_df <- proxy_joined %>%
    select(core_base, y_bp, all_of(resp), all_of(top_preds)) %>%
    pivot_longer(cols = all_of(top_preds), names_to = "predictor", values_to = "x") %>%
    mutate(
      y = as_numeric_safe(.data[[resp]]),
      x = as_numeric_safe(x)
    ) %>%
    filter(!is.na(x), !is.na(y))

  if (nrow(plot_df) == 0) return(NULL)

  p <- ggplot(plot_df, aes(x = x, y = y, color = core_base)) +
    geom_point(alpha = 0.7, size = 1.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    facet_wrap(~ predictor, scales = "free_x") +
    labs(
      title = paste0("Top proxy associations for ", resp),
      subtitle = "Points colored by core; black line = per-panel linear fit",
      x = "Proxy value",
      y = resp,
      color = "Core"
    ) +
    theme_bw(base_size = 10)

  out_plot <- file.path(plots_dir, paste0("top_associations_", resp, ".png"))
  ggsave(out_plot, p, width = 13, height = 8, dpi = 220)
})

summary_table <- association_results %>%
  group_by(response) %>%
  summarise(
    n_predictors_tested = sum(!is.na(spearman_rho)),
    n_sig_spearman_fdr05 = sum(!is.na(spearman_fdr) & spearman_fdr <= 0.05),
    n_sig_lm_fdr05 = sum(!is.na(lm_fdr) & lm_fdr <= 0.05),
    top_predictor_by_abs_rho = predictor[which.max(abs(replace_na(spearman_rho, 0)))],
    top_abs_spearman_rho = max(abs(spearman_rho), na.rm = TRUE),
    .groups = "drop"
  )

metadata_out <- tibble(
  model_version = params$model_version,
  run_timestamp = as.character(Sys.time()),
  sample_baselines_file = input_sample_baselines,
  metadata_file = input_metadata,
  foraminifera_file = input_forams,
  sst_file = input_sst,
  xrf_file = input_xrf,
  geob25202_file = input_geob,
  depth_match_tolerance_cm = params$depth_match_tolerance_cm,
  min_pairs_for_test = params$min_pairs_for_test,
  collapsed_rows = nrow(damage_agg),
  predictors_considered = length(predictor_vars),
  random_seed = params$random_seed,
  assumption_core_mapping = "core_base computed by removing _R1/_R2 from metadata core",
  assumption_depth_matching = "Nearest depth match within tolerance by core_base"
)

proxy_joined_file <- file.path(out_dir, "damage-baselines-with-proxies.tsv.gz")
associations_file <- file.path(out_dir, "proxy-associations.tsv")
top_hits_file <- file.path(out_dir, "proxy-associations-top-hits.tsv")
summary_file <- file.path(out_dir, "summary.tsv")
metadata_file <- file.path(out_dir, "metadata.tsv")
session_file <- file.path(out_dir, "sessioninfo.txt")

log_msg("Writing outputs")
write_tsv(proxy_joined, proxy_joined_file)
write_tsv(association_results, associations_file)
write_tsv(top_hits, top_hits_file)
write_tsv(summary_table, summary_file)
write_tsv(metadata_out, metadata_file)
writeLines(capture.output(sessionInfo()), session_file)

checksum_targets <- c(
  proxy_joined_file,
  associations_file,
  top_hits_file,
  summary_file,
  metadata_file,
  session_file
)

plot_files <- list.files(plots_dir, pattern = "\\.png$", full.names = TRUE)
checksum_targets <- c(checksum_targets, plot_files)

checksums <- tibble(
  file = basename(checksum_targets),
  md5 = unname(tools::md5sum(checksum_targets))
)
checksums_file <- file.path(out_dir, "output-checksums.tsv")
write_tsv(checksums, checksums_file)

log_msg("Completed successfully")
