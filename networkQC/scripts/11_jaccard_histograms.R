#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

OUT_FULL <- here("networkQC", "results", "full_eval")
OUT_FIG <- here("networkQC", "results", "figures")
OUT_TAB <- here("networkQC", "results", "tables")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

candidate_settings <- c("baseline", "exp1", "exp2", "exp3", "exp4", "opt1", "opt2", "opt3", "opt4", "opt5")

read_bootstrap <- function(setting_id) {
  path <- file.path(OUT_FULL, setting_id, "bootstrap_module_stability.tsv")
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (!"jaccard" %in% names(dt) || !"module" %in% names(dt)) return(NULL)
  dt[, setting_id := setting_id]
  dt[is.finite(jaccard)]
}

boot_dt <- rbindlist(lapply(candidate_settings, read_bootstrap), fill = TRUE)
if (!nrow(boot_dt)) {
  stop("No bootstrap Jaccard files found for requested settings.")
}

present_settings <- unique(boot_dt$setting_id)
missing_settings <- setdiff(candidate_settings, present_settings)
if (length(missing_settings)) {
  message("Missing settings (skipped): ", paste(missing_settings, collapse = ", "))
}

boot_dt[, setting_id := factor(setting_id, levels = candidate_settings[candidate_settings %in% present_settings])]

sum_dt <- boot_dt[, .(
  n = .N,
  jaccard_mean = mean(jaccard, na.rm = TRUE),
  jaccard_median = median(jaccard, na.rm = TRUE),
  jaccard_p05 = as.numeric(quantile(jaccard, 0.05, na.rm = TRUE)),
  jaccard_p95 = as.numeric(quantile(jaccard, 0.95, na.rm = TRUE))
), by = .(setting_id, module)]
setorder(sum_dt, setting_id, -jaccard_median)
fwrite(sum_dt, file.path(OUT_TAB, "jaccard_histogram_summary_by_setting_module.tsv"), sep = "\t")

overall_dt <- boot_dt[, .(
  n = .N,
  jaccard_mean = mean(jaccard, na.rm = TRUE),
  jaccard_median = median(jaccard, na.rm = TRUE),
  jaccard_p05 = as.numeric(quantile(jaccard, 0.05, na.rm = TRUE)),
  jaccard_p95 = as.numeric(quantile(jaccard, 0.95, na.rm = TRUE))
), by = setting_id]
setorder(overall_dt, setting_id)
fwrite(overall_dt, file.path(OUT_TAB, "jaccard_histogram_summary_by_setting.tsv"), sep = "\t")

exp3_dt <- boot_dt[setting_id == "exp3"]
if (nrow(exp3_dt)) {
  p_exp3_all <- ggplot(exp3_dt, aes(x = jaccard)) +
    geom_histogram(binwidth = 0.025, boundary = 0, fill = "#1f77b4", color = "white", alpha = 0.95) +
    geom_vline(xintercept = mean(exp3_dt$jaccard, na.rm = TRUE), color = "#d62728", linewidth = 0.8) +
    geom_vline(xintercept = median(exp3_dt$jaccard, na.rm = TRUE), color = "#2ca02c", linewidth = 0.8, linetype = "dashed") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    labs(
      title = "exp3 Bootstrap Jaccard Distribution",
      subtitle = sprintf(
        "n=%d | mean=%.3f | median=%.3f | p05=%.3f | p95=%.3f",
        nrow(exp3_dt),
        mean(exp3_dt$jaccard, na.rm = TRUE),
        median(exp3_dt$jaccard, na.rm = TRUE),
        as.numeric(quantile(exp3_dt$jaccard, 0.05, na.rm = TRUE)),
        as.numeric(quantile(exp3_dt$jaccard, 0.95, na.rm = TRUE))
      ),
      x = "Best-match Jaccard overlap",
      y = "Bootstrap module matches"
    ) +
    theme_minimal(base_size = 11)

  p_exp3_by_module <- ggplot(exp3_dt, aes(x = jaccard)) +
    geom_histogram(binwidth = 0.025, boundary = 0, fill = "#4c78a8", color = "white", alpha = 0.95) +
    facet_wrap(~ module, ncol = 4, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = "exp3 Bootstrap Jaccard by Module",
      subtitle = "Each panel is one reference module; values are best-match overlaps across bootstrap reruns",
      x = "Jaccard overlap",
      y = "Count"
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))

  ggsave(
    filename = file.path(OUT_FIG, "jaccard_hist_exp3_overall.png"),
    plot = p_exp3_all,
    width = 8.8,
    height = 5.6,
    dpi = 220
  )
  ggsave(
    filename = file.path(OUT_FIG, "jaccard_hist_exp3_by_module.png"),
    plot = p_exp3_by_module,
    width = 12.5,
    height = 7.8,
    dpi = 220
  )
}

order_dt <- fread(file.path(OUT_FULL, "all_settings_ranked.tsv"))
order_vec <- order_dt$setting_id[order_dt$setting_id %in% levels(boot_dt$setting_id)]
if (length(order_vec)) {
  boot_dt[, setting_id := factor(as.character(setting_id), levels = order_vec)]
}

p_all <- ggplot(boot_dt, aes(x = jaccard)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "#6baed6", color = "white", alpha = 0.95) +
  facet_wrap(~ setting_id, ncol = 4, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Bootstrap Jaccard Distributions Across Candidate Settings",
    subtitle = "Settings: baseline, exp1-4, opt1-5",
    x = "Best-match Jaccard overlap",
    y = "Count"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

ggsave(
  filename = file.path(OUT_FIG, "jaccard_hist_candidates_panel.png"),
  plot = p_all,
  width = 14.8,
  height = 9.6,
  dpi = 220
)

message("[networkQC] wrote Jaccard histogram figures and summary tables")
