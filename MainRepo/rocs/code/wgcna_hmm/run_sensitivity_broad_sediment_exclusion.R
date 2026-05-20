#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))

ensure_dirs(c(DIRS$results, DIRS$main, DIRS$reports, DIRS$logs))

out_dir <- file.path(DIRS$results, "sensitivity_broad_sediment_exclusion")
baseline_dir <- file.path(out_dir, "baseline")
sensitivity_dir <- file.path(out_dir, "sensitivity")
ensure_dirs(c(out_dir, baseline_dir, sensitivity_dir))

extra_phyla <- c("p__Aerophobota", "p__Atribacterota", "p__Asgardarchaeota")
extra_phyla_csv <- paste(extra_phyla, collapse = ",")

key_files <- c(
  key_metrics = file.path("results", "wgcna_hmm", "reports", "key_metrics.tsv"),
  taxa_after_filter = file.path("results", "wgcna_hmm", "main", "taxa_after_filter.tsv"),
  module_assignments = file.path("results", "wgcna_hmm", "main", "module_assignments.tsv"),
  module_sizes = file.path("results", "wgcna_hmm", "main", "module_sizes.tsv"),
  hmm_states_main = file.path("results", "wgcna_hmm", "main", "hmm_states_main.tsv"),
  state_concordance = file.path("results", "wgcna_hmm", "main", "state_concordance_R1_vs_R2_summary.tsv"),
  module_xrf_results = file.path("results", "wgcna_hmm", "main", "module_xrf_results.tsv"),
  loco_module_xrf = file.path("results", "wgcna_hmm", "main", "loco_module_xrf.tsv"),
  tea_module_enrichment = file.path("results", "wgcna_hmm", "main", "tea", "tea_module_enrichment_centered_wilcoxon.tsv"),
  tea_state_summary = file.path("results", "wgcna_hmm", "main", "tea", "tea_state_summary_main.tsv")
)

sensitivity_extra_files <- c(
  data_prep_summary = file.path("results", "wgcna_hmm", "main", "data_prep_summary.tsv"),
  data_prep_run_metadata = file.path("results", "wgcna_hmm", "main", "01_data_prep_run_metadata.tsv"),
  deep_manifest = file.path("results", "wgcna_hmm", "main", "deep_sediment_excluded_taxa.tsv")
)

copy_named_files <- function(named_rel_paths, dest_dir) {
  for (nm in names(named_rel_paths)) {
    src <- file.path(REPO_ROOT, named_rel_paths[[nm]])
    if (!file.exists(src)) {
      stop("Required file missing: ", src)
    }
    dest <- file.path(dest_dir, basename(named_rel_paths[[nm]]))
    ok <- file.copy(src, dest, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to copy file: ", src, " -> ", dest)
    }
  }
}

run_workflow <- function(env = character()) {
  cmd <- file.path("code", "wgcna_hmm", "run_workflow.R")
  status <- system2("Rscript", args = cmd, env = env)
  if (!identical(status, 0L)) {
    stop("Workflow execution failed with status ", status)
  }
}

pairwise_partition_agreement_pct <- function(a, b) {
  n <- length(a)
  if (n <= 1) return(100)
  idx <- utils::combn(seq_len(n), 2)
  same_a <- a[idx[1, ]] == a[idx[2, ]]
  same_b <- b[idx[1, ]] == b[idx[2, ]]
  mean(same_a == same_b) * 100
}

log_msg("Step 1/5: Snapshotting current baseline outputs...")
copy_named_files(key_files, baseline_dir)

log_msg("Step 2/5: Running broadened exclusion sensitivity workflow...")
run_workflow(env = c(
  "WGCNA_HMM_EXCLUDE_DEEP_SEDIMENT_TAXA=1",
  paste0("WGCNA_HMM_EXTRA_EXCLUDED_PHYLA=", extra_phyla_csv)
))

log_msg("Step 3/5: Capturing sensitivity artifacts...")
copy_named_files(key_files, sensitivity_dir)
copy_named_files(sensitivity_extra_files, sensitivity_dir)

log_msg("Step 4/5: Building comparison artifacts...")

baseline_key <- fread(file.path(baseline_dir, "key_metrics.tsv"))
sens_key <- fread(file.path(sensitivity_dir, "key_metrics.tsv"))

metrics_keep <- c(
  "n_taxa_main", "n_modules_non_grey", "hmm_states_main", "r1_r2_state_pct_match",
  "xrf_significant_n", "loco_stable_n", "tea_module_enrichment_calls_n"
)

comparison_key <- merge(
  baseline_key[metric %in% metrics_keep, .(metric, value_baseline = as.numeric(value))],
  sens_key[metric %in% metrics_keep, .(metric, value_sensitivity = as.numeric(value))],
  by = "metric",
  all = TRUE
)
comparison_key[, delta_sensitivity_minus_baseline := value_sensitivity - value_baseline]
fwrite(comparison_key, file.path(out_dir, "comparison_key_metrics.tsv"), sep = "\t")

manifest <- fread(file.path(sensitivity_dir, "deep_sediment_excluded_taxa.tsv"))
dprep <- fread(file.path(sensitivity_dir, "data_prep_summary.tsv"))

by_reason <- manifest[, .(excluded_taxa_n = uniqueN(taxon)), by = exclusion_reason][order(-excluded_taxa_n, exclusion_reason)]
fwrite(by_reason, file.path(out_dir, "comparison_exclusion_by_reason.tsv"), sep = "\t")

by_reason_phylum <- manifest[, .(excluded_taxa_n = uniqueN(taxon)), by = .(exclusion_reason, phylum)][order(exclusion_reason, -excluded_taxa_n, phylum)]
fwrite(by_reason_phylum, file.path(out_dir, "comparison_exclusion_by_reason_phylum.tsv"), sep = "\t")

summary_rows <- rbindlist(list(
  data.table(metric = "excluded_taxa_n_manifest", value = as.character(uniqueN(manifest$taxon))),
  data.table(metric = "excluded_manifest_rows", value = as.character(nrow(manifest))),
  data.table(metric = "excluded_taxa_n_data_prep_summary", value = dprep[metric == "deep_sediment_excluded_selected_samples", value][1]),
  data.table(metric = "excluded_taxa_prevalent_main_window_n", value = dprep[metric == "deep_sediment_excluded_prevalent_main_window", value][1]),
  data.table(metric = "deep_sediment_filter_enabled", value = dprep[metric == "deep_sediment_filter_enabled", value][1]),
  data.table(metric = "extra_excluded_phyla_count", value = dprep[metric == "extra_excluded_phyla_count", value][1]),
  data.table(metric = "extra_excluded_phyla_values", value = dprep[metric == "extra_excluded_phyla_values", value][1]),
  data.table(metric = "extra_phylum_candidates_total", value = dprep[metric == "extra_phylum_candidates_total", value][1]),
  data.table(metric = "extra_phylum_excluded_selected_samples", value = dprep[metric == "extra_phylum_excluded_selected_samples", value][1]),
  data.table(metric = "extra_phylum_excluded_prevalent_main_window", value = dprep[metric == "extra_phylum_excluded_prevalent_main_window", value][1])
))
fwrite(summary_rows, file.path(out_dir, "comparison_exclusion_summary.tsv"), sep = "\t")

base_mod <- fread(file.path(baseline_dir, "module_assignments.tsv"))
sens_mod <- fread(file.path(sensitivity_dir, "module_assignments.tsv"))
mods <- merge(
  base_mod[, .(taxon, module_baseline = module)],
  sens_mod[, .(taxon, module_sensitivity = module)],
  by = "taxon",
  all = FALSE
)
mods[, changed := module_baseline != module_sensitivity]
mod_summary <- data.table(
  metric = c("shared_taxa_compared", "changed_module_assignment_n", "changed_module_assignment_pct"),
  value = c(nrow(mods), sum(mods$changed), if (nrow(mods) > 0) 100 * mean(mods$changed) else 0)
)
fwrite(mod_summary, file.path(out_dir, "comparison_module_membership_summary.tsv"), sep = "\t")
mod_cont <- mods[, .N, by = .(module_baseline, module_sensitivity)][order(module_baseline, module_sensitivity)]
fwrite(mod_cont, file.path(out_dir, "comparison_module_membership_contingency.tsv"), sep = "\t")

base_hmm <- fread(file.path(baseline_dir, "hmm_states_main.tsv"))
sens_hmm <- fread(file.path(sensitivity_dir, "hmm_states_main.tsv"))
hmm_cmp <- merge(
  base_hmm[, .(sample, state_baseline = state, label_baseline = label)],
  sens_hmm[, .(sample, state_sensitivity = state, label_sensitivity = label)],
  by = "sample",
  all = FALSE
)
hmm_cmp[, exact_state_label_match := state_baseline == state_sensitivity]
fwrite(hmm_cmp, file.path(out_dir, "comparison_hmm_states_by_sample.tsv"), sep = "\t")

hmm_summary <- data.table(
  metric = c("samples_compared", "exact_state_label_match_pct", "pairwise_partition_agreement_pct", "note"),
  value = c(
    nrow(hmm_cmp),
    if (nrow(hmm_cmp) > 0) 100 * mean(hmm_cmp$exact_state_label_match) else 0,
    pairwise_partition_agreement_pct(hmm_cmp$state_baseline, hmm_cmp$state_sensitivity),
    "State IDs may be permuted across runs; partition agreement is label-invariant."
  )
)
fwrite(hmm_summary, file.path(out_dir, "comparison_hmm_state_summary.tsv"), sep = "\t")

base_tea <- fread(file.path(baseline_dir, "tea_module_enrichment_centered_wilcoxon.tsv"))
sens_tea <- fread(file.path(sensitivity_dir, "tea_module_enrichment_centered_wilcoxon.tsv"))
tea_cmp <- merge(
  base_tea[, .(module, index, analysis_set, direction_baseline = direction, q_value_baseline = q_value, status_baseline = status)],
  sens_tea[, .(module, index, analysis_set, direction_sensitivity = direction, q_value_sensitivity = q_value, status_sensitivity = status)],
  by = c("module", "index", "analysis_set"),
  all = TRUE
)
fwrite(tea_cmp, file.path(out_dir, "comparison_tea_module_enrichment.tsv"), sep = "\t")

brown_check <- tea_cmp[module == "brown" & index == "OAP" & analysis_set == "include_grey"]
brown_line <- if (nrow(brown_check) == 1) {
  sprintf(
    "brown OAP include_grey: baseline=%s (q=%s), sensitivity=%s (q=%s)",
    brown_check$direction_baseline,
    format(brown_check$q_value_baseline, scientific = TRUE),
    brown_check$direction_sensitivity,
    format(brown_check$q_value_sensitivity, scientific = TRUE)
  )
} else {
  "brown OAP include_grey row missing in one or both runs"
}

summary_lines <- c(
  "# Broad deep-sediment exclusion sensitivity comparison",
  "",
  "## Exclusion definition used",
  "- curated diagenetic rule:",
  "  - `display_category == \"Diagenetic\"`",
  "  - `signal_source %in% c(\"diagenetic_in_situ\", \"benthic_resident\")`",
  "  - `confidence_score >= 0.85`",
  sprintf("- plus extra excluded phyla: `%s`", paste(extra_phyla, collapse = "`, `")),
  "",
  "## Exclusion counts",
  sprintf("- Excluded taxa (unique): %s", summary_rows[metric == "excluded_taxa_n_manifest", value][1]),
  sprintf("- Excluded manifest rows (taxon x reason): %s", summary_rows[metric == "excluded_manifest_rows", value][1]),
  sprintf("- Excluded taxa that would pass prevalence in main window: %s", summary_rows[metric == "excluded_taxa_prevalent_main_window_n", value][1]),
  "",
  "## Key metric deltas (sensitivity - baseline)",
  paste(sprintf(
    "- %s: baseline=%s, sensitivity=%s, delta=%s",
    comparison_key$metric,
    comparison_key$value_baseline,
    comparison_key$value_sensitivity,
    comparison_key$delta_sensitivity_minus_baseline
  ), collapse = "\n"),
  "",
  "## Module membership stability",
  sprintf("- Shared taxa compared: %s", mod_summary[metric == "shared_taxa_compared", value][1]),
  sprintf("- Taxa changing assigned module: %s (%.2f%%)", mod_summary[metric == "changed_module_assignment_n", value][1], mod_summary[metric == "changed_module_assignment_pct", value][1]),
  "",
  "## HMM state comparison",
  sprintf("- Exact state label match across samples: %.2f%%", as.numeric(hmm_summary[metric == "exact_state_label_match_pct", value][1])),
  sprintf("- Label-invariant pairwise partition agreement: %.2f%%", as.numeric(hmm_summary[metric == "pairwise_partition_agreement_pct", value][1])),
  "",
  "## Brown-associated TEA pattern check",
  paste0("- ", brown_line),
  "",
  "## Interpretation (conservative)",
  "- Architecture-level metrics are compared descriptively; directional TEA changes indicate sensitivity for specific module-index associations.",
  "- These diagnostics are not causal attribution and should be interpreted with classification uncertainty in mind."
 )
writeLines(summary_lines, con = file.path(out_dir, "comparison_summary.md"))

log_msg("Step 5/5: Restoring default outputs by rerunning workflow with defaults...")
run_workflow(env = c(
  "WGCNA_HMM_EXCLUDE_DEEP_SEDIMENT_TAXA=0",
  "WGCNA_HMM_EXTRA_EXCLUDED_PHYLA="
))

restore_md <- fread(file.path(DIRS$main, "01_data_prep_run_metadata.tsv"))
restore_checks <- data.table(
  check = c(
    "exclude_deep_sediment_taxa_default_off",
    "extra_excluded_phyla_default_empty"
  ),
  value = c(
    restore_md$exclude_deep_sediment_taxa[1],
    restore_md$extra_excluded_phyla[1]
  )
)
fwrite(restore_checks, file.path(out_dir, "restore_defaults_check.tsv"), sep = "\t")

log_msg("Broad sensitivity run complete. Outputs: ", out_dir)
