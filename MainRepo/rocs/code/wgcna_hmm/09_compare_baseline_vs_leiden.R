#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))

baseline_root <- file.path(REPO_ROOT, "results", "wgcna_hmm")
alt_root <- DIRS$results

if (normalizePath(baseline_root, winslash = "/", mustWork = FALSE) == normalizePath(alt_root, winslash = "/", mustWork = FALSE)) {
  stop("Comparison script requires an alternative results root distinct from baseline results/wgcna_hmm.")
}

cmp_dir <- file.path(alt_root, "comparison_vs_baseline")
ensure_dirs(c(cmp_dir))

read_key_metrics <- function(root) {
  path <- file.path(root, "reports", "key_metrics.tsv")
  if (!file.exists(path)) stop("Missing key metrics: ", path)
  fread(path)
}

read_safe <- function(path) {
  if (!file.exists(path)) return(data.table())
  fread(path)
}

base_key <- read_key_metrics(baseline_root)
alt_key <- read_key_metrics(alt_root)

metrics_keep <- c(
  "n_taxa_main",
  "n_modules_non_grey",
  "hmm_states_main",
  "r1_r2_state_pct_match",
  "xrf_significant_n",
  "loco_stable_n",
  "tea_module_enrichment_calls_n"
)

key_cmp <- merge(
  base_key[metric %in% metrics_keep, .(metric, baseline = as.character(value))],
  alt_key[metric %in% metrics_keep, .(metric, leiden = as.character(value))],
  by = "metric",
  all = TRUE
)

key_cmp[, `:=`(
  baseline_num = suppressWarnings(as.numeric(baseline)),
  leiden_num = suppressWarnings(as.numeric(leiden))
)]
key_cmp[, delta_leiden_minus_baseline := leiden_num - baseline_num]
fwrite(key_cmp, file.path(cmp_dir, "comparison_key_metrics.tsv"), sep = "\t")

base_mod <- read_safe(file.path(baseline_root, "main", "module_assignments.tsv"))
alt_mod <- read_safe(file.path(alt_root, "main", "module_assignments.tsv"))
base_sizes <- read_safe(file.path(baseline_root, "main", "module_sizes.tsv"))
alt_sizes <- read_safe(file.path(alt_root, "main", "module_sizes.tsv"))

if (nrow(base_mod) == 0 || nrow(alt_mod) == 0) {
  stop("Missing module assignment files in baseline or Leiden outputs.")
}

grey_summary <- data.table(
  run = c("baseline", "leiden"),
  n_taxa_main = c(nrow(base_mod), nrow(alt_mod)),
  grey_size = c(sum(base_mod$module == "grey"), sum(alt_mod$module == "grey"))
)
grey_summary[, grey_fraction := grey_size / n_taxa_main]
fwrite(grey_summary, file.path(cmp_dir, "comparison_grey_burden.tsv"), sep = "\t")

summarize_sizes <- function(sizes_dt, run_name) {
  if (nrow(sizes_dt) == 0) return(data.table())
  ng <- sizes_dt[module != "grey", n_taxa]
  data.table(
    run = run_name,
    n_modules_non_grey = sum(sizes_dt$module != "grey"),
    min_non_grey = if (length(ng) > 0) min(ng) else NA_real_,
    q25_non_grey = if (length(ng) > 0) stats::quantile(ng, probs = 0.25, names = FALSE) else NA_real_,
    median_non_grey = if (length(ng) > 0) stats::median(ng) else NA_real_,
    q75_non_grey = if (length(ng) > 0) stats::quantile(ng, probs = 0.75, names = FALSE) else NA_real_,
    max_non_grey = if (length(ng) > 0) max(ng) else NA_real_
  )
}

size_dist <- rbindlist(list(
  summarize_sizes(base_sizes, "baseline"),
  summarize_sizes(alt_sizes, "leiden")
), use.names = TRUE, fill = TRUE)
fwrite(size_dist, file.path(cmp_dir, "comparison_module_size_distribution_summary.tsv"), sep = "\t")

base_pres <- read_safe(file.path(baseline_root, "main", "module_preservation_validation.tsv"))
alt_pres <- read_safe(file.path(alt_root, "main", "module_preservation_validation.tsv"))
pres_counts <- rbindlist(list(
  if (nrow(base_pres) > 0) base_pres[, .N, by = preserved][, run := "baseline"] else data.table(),
  if (nrow(alt_pres) > 0) alt_pres[, .N, by = preserved][, run := "leiden"] else data.table()
), fill = TRUE)
if (nrow(pres_counts) > 0) fwrite(pres_counts[order(run, preserved)], file.path(cmp_dir, "comparison_module_preservation_counts.tsv"), sep = "\t")

summary_lines <- c(
  "# Baseline vs Leiden alternative workflow comparison",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Leiden modules here are a graph-partitioning sensitivity analysis and are not directly equivalent",
  "to dynamic-tree modules from consensus WGCNA.",
  "",
  "## Grey burden",
  sprintf("- Baseline grey: %d / %d (%.2f%%)",
          grey_summary[run == "baseline", grey_size],
          grey_summary[run == "baseline", n_taxa_main],
          100 * grey_summary[run == "baseline", grey_fraction]),
  sprintf("- Leiden grey: %d / %d (%.2f%%)",
          grey_summary[run == "leiden", grey_size],
          grey_summary[run == "leiden", n_taxa_main],
          100 * grey_summary[run == "leiden", grey_fraction]),
  sprintf("- Delta grey fraction (Leiden - baseline): %.2f percentage points",
          100 * (grey_summary[run == "leiden", grey_fraction] - grey_summary[run == "baseline", grey_fraction])),
  "",
  "## Key downstream metrics",
  paste(sprintf("- %s: baseline=%s; leiden=%s; delta=%s",
                key_cmp$metric,
                key_cmp$baseline,
                key_cmp$leiden,
                ifelse(is.finite(key_cmp$delta_leiden_minus_baseline),
                       as.character(key_cmp$delta_leiden_minus_baseline), "NA")),
        collapse = "\n"),
  "",
  "## Interpretation",
  "- Check whether grey burden is reduced while preserving high-level downstream signals.",
  "- Strong shifts in HMM/XRF/TEA metrics indicate sensitivity to module-construction choice."
)

writeLines(summary_lines, con = file.path(cmp_dir, "comparison_summary.md"))

write_run_metadata(
  file.path(cmp_dir, "09_compare_baseline_vs_leiden_run_metadata.tsv"),
  "09_compare_baseline_vs_leiden.R",
  extra = list(
    baseline_root = baseline_root,
    leiden_root = alt_root,
    module_method = PARAMS$module_method,
    leiden_resolution = PARAMS$leiden_resolution
  )
)
write_session_info(file.path(cmp_dir, "09_compare_baseline_vs_leiden_sessionInfo.txt"))

log_msg("09_compare_baseline_vs_leiden complete.")
