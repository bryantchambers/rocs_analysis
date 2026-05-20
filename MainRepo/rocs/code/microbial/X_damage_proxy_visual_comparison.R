library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(tibble)

ROCS_ROOT <- Sys.getenv("ROCS_ROOT", "/projects/caeg/people/ngm902/apps/repos/rocs")
ROCS_ROOT <- normalizePath(ROCS_ROOT, mustWork = FALSE)

params <- list(
  model_version = "damage-proxy-visual-comparison-v1",
  random_seed = 20260325,
  selected_vars = c(
    "sample_A_b_median",
    "xrf_nd",
    "xrf_mn",
    "xrf_ca_fe",
    "xrf_rb_sr",
    "xrf_ti_ca",
    "sst_sst_uk37_alkenone",
    "foram_g_bulloides_pct"
  )
)

set.seed(params$random_seed)

input_file <- file.path(
  ROCS_ROOT,
  "results/microbial/damage/damage-baseline-proxy-associations/damage-baselines-with-proxies.tsv.gz"
)

out_dir <- file.path(
  ROCS_ROOT,
  "results/microbial/damage/damage-baseline-proxy-associations/plots/visual-comparison"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")

log_msg <- function(...) {
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = log_file, append = TRUE)
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - m) / s
}

pretty_labels <- c(
  sample_A_b_median = "Damage baseline (sample_A_b_median)",
  xrf_nd = "XRF Nd",
  xrf_mn = "XRF Mn",
  xrf_ca_fe = "XRF Ca/Fe",
  xrf_rb_sr = "XRF Rb/Sr",
  xrf_ti_ca = "XRF Ti/Ca",
  sst_sst_uk37_alkenone = "SST UK'37 alkenone",
  foram_g_bulloides_pct = "Foram G. bulloides (%)"
)

log_msg("Reading merged damage-proxy table: ", input_file)
dat <- read_tsv(input_file, show_col_types = FALSE)

required_cols <- c("core_base", "y_bp")
missing_required <- setdiff(required_cols, names(dat))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

present_vars <- params$selected_vars[params$selected_vars %in% names(dat)]
missing_vars <- setdiff(params$selected_vars, present_vars)
if (length(missing_vars) > 0) {
  log_msg("Warning: missing selected variables and skipped: ", paste(missing_vars, collapse = ", "))
}

plot_long <- dat %>%
  select(core_base, y_bp, all_of(present_vars)) %>%
  mutate(y_bp = as.numeric(y_bp)) %>%
  pivot_longer(cols = all_of(present_vars), names_to = "variable", values_to = "value") %>%
  filter(!is.na(core_base), !is.na(y_bp), is.finite(y_bp)) %>%
  group_by(core_base, variable) %>%
  mutate(
    value_z = zscore(as.numeric(value)),
    n_non_na = sum(!is.na(value))
  ) %>%
  ungroup() %>%
  mutate(variable_label = pretty_labels[variable])

plot_long <- plot_long %>%
  filter(!is.na(variable_label), !is.na(value_z), is.finite(value_z), n_non_na >= 5)

if (nrow(plot_long) == 0) {
  stop("No plottable rows after filtering. Check selected variables and data coverage.")
}

log_msg("Generating panel plot by core and variable")
p_panel <- ggplot(plot_long, aes(x = y_bp, y = value_z, group = interaction(core_base, variable))) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_line(color = "#2c3e50", linewidth = 0.35, alpha = 0.9) +
  geom_point(color = "#2c3e50", size = 0.55, alpha = 0.75) +
  facet_grid(variable_label ~ core_base, scales = "free_x") +
  scale_x_reverse() +
  labs(
    title = "Visual comparison of damage baseline and selected proxies",
    subtitle = "Each panel is z-scored within core and variable (mean=0, sd=1)",
    x = "Age (y BP)",
    y = "Standardized value (z-score)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0),
    plot.title = element_text(face = "bold")
  )

panel_file <- file.path(out_dir, "damage_proxy_panel_zscore_by_core.png")
ggsave(panel_file, p_panel, width = 15, height = 12, dpi = 260)

log_msg("Generating overlay plot for fast cross-core comparison")
overlay_vars <- c("sample_A_b_median", "xrf_nd", "xrf_mn", "foram_g_bulloides_pct", "sst_sst_uk37_alkenone")
overlay_vars <- overlay_vars[overlay_vars %in% present_vars]

plot_overlay <- plot_long %>%
  filter(variable %in% overlay_vars)

p_overlay <- ggplot(plot_overlay, aes(x = y_bp, y = value_z, color = core_base)) +
  geom_line(linewidth = 0.35, alpha = 0.8) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.65) +
  facet_wrap(~ variable_label, ncol = 1, scales = "free_x") +
  scale_x_reverse() +
  labs(
    title = "Damage and proxy trends by core",
    subtitle = "Thin lines = observed; thick lines = loess trend (z-score scale)",
    x = "Age (y BP)",
    y = "Standardized value (z-score)",
    color = "Core"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

overlay_file <- file.path(out_dir, "damage_proxy_overlay_trends.png")
ggsave(overlay_file, p_overlay, width = 11, height = 13, dpi = 260)

metadata <- tibble(
  model_version = params$model_version,
  run_timestamp = as.character(Sys.time()),
  input_file = input_file,
  random_seed = params$random_seed,
  selected_vars_requested = paste(params$selected_vars, collapse = ";"),
  selected_vars_present = paste(present_vars, collapse = ";"),
  panel_plot_file = basename(panel_file),
  overlay_plot_file = basename(overlay_file),
  zscore_definition = "Within each core_base and variable"
)

metadata_file <- file.path(out_dir, "metadata.tsv")
session_file <- file.path(out_dir, "sessioninfo.txt")
checksums_file <- file.path(out_dir, "output-checksums.tsv")

write_tsv(metadata, metadata_file)
writeLines(capture.output(sessionInfo()), session_file)

checksum_targets <- c(panel_file, overlay_file, metadata_file, session_file)
checksums <- tibble(
  file = basename(checksum_targets),
  md5 = unname(tools::md5sum(checksum_targets))
)
write_tsv(checksums, checksums_file)

log_msg("Completed successfully")
