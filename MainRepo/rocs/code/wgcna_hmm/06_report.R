#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))

meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
taxa <- fread(file.path(DIRS$main, "taxa_after_filter.tsv"))
mods <- fread(file.path(DIRS$main, "module_assignments.tsv"))
pres <- fread(file.path(DIRS$main, "module_preservation_validation.tsv"))
hmm <- fread(file.path(DIRS$main, "hmm_states_main.tsv"))
balance_summary <- if (file.exists(file.path(DIRS$main, "balance_design_summary.tsv"))) {
  fread(file.path(DIRS$main, "balance_design_summary.tsv"))
} else {
  data.table(metric = character(), value = character())
}
state_fp <- if (file.exists(file.path(DIRS$main, "state_fingerprints_main.tsv"))) {
  fread(file.path(DIRS$main, "state_fingerprints_main.tsv"))
} else {
  data.table()
}
concord <- fread(file.path(DIRS$main, "state_concordance_R1_vs_R2_summary.tsv"))
xrf_res <- fread(file.path(DIRS$main, "module_xrf_results.tsv"))
loco <- fread(file.path(DIRS$main, "loco_module_xrf.tsv"))
mode_eval <- fread(file.path(DIRS$extended, "extended_harmonization_mode_evaluation.tsv"))
ext_diag <- fread(file.path(DIRS$extended, "extended_older_period_diagnostics.tsv"))
ext_states <- fread(file.path(DIRS$extended, "hmm_states_extended_selected.tsv"))
ext_scale_diag <- fread(file.path(DIRS$extended, "extended_corewise_z_scaling_fallback_diagnostics.tsv"))
ext_state_origin <- fread(file.path(DIRS$extended, "extended_state_origin_summary.tsv"))
ext_lock_diag <- fread(file.path(DIRS$extended, "extended_locked_overlap_diagnostics.tsv"))

tea_dir <- file.path(DIRS$main, "tea")
tea_sample <- if (file.exists(file.path(tea_dir, "tea_indices_per_sample_main.tsv"))) {
  fread(file.path(tea_dir, "tea_indices_per_sample_main.tsv"))
} else {
  data.table()
}
tea_mod_enrich <- if (file.exists(file.path(tea_dir, "tea_module_enrichment_centered_wilcoxon.tsv"))) {
  fread(file.path(tea_dir, "tea_module_enrichment_centered_wilcoxon.tsv"))
} else {
  data.table()
}
tea_taxon <- if (file.exists(file.path(tea_dir, "tea_taxon_metabolic_calls.tsv"))) {
  fread(file.path(tea_dir, "tea_taxon_metabolic_calls.tsv"))
} else {
  data.table()
}
tea_taxa_with_oap_n <- if ("has_oap" %in% names(tea_taxon)) nrow(tea_taxon[has_oap == TRUE]) else 0L
get_balance_value <- function(metric_name, default = NA_character_) {
  val <- balance_summary[metric == metric_name, value]
  if (length(val) == 0) return(default)
  as.character(val[1])
}

nsamp <- nrow(meta_main)
ntaxa <- nrow(taxa)
nmods <- uniqueN(mods[module != "grey", module])

taxa_per_module <- mods[, .(n_taxa = .N), by = module][order(-n_taxa)]
fwrite(taxa_per_module, file.path(DIRS$reports, "taxa_per_module.tsv"), sep = "\t")

module_taxa <- mods[order(module, taxon)]
fwrite(module_taxa, file.path(DIRS$reports, "module_taxa_membership.tsv"), sep = "\t")

pres_summary <- pres[, .(
  n_modules = .N,
  strong = sum(preserved == "strong", na.rm = TRUE),
  moderate = sum(preserved == "moderate", na.rm = TRUE),
  weak = sum(preserved == "weak", na.rm = TRUE),
  median_Zsummary = median(Zsummary, na.rm = TRUE)
)]
fwrite(pres_summary, file.path(DIRS$reports, "module_preservation_summary.tsv"), sep = "\t")

hmm_summary <- hmm[, .(
  n_samples = .N,
  n_states = uniqueN(state),
  min_age_kyr = min(age_kyr, na.rm = TRUE),
  max_age_kyr = max(age_kyr, na.rm = TRUE)
), by = core]
fwrite(hmm_summary, file.path(DIRS$reports, "hmm_state_summary_by_core.tsv"), sep = "\t")

if (nrow(state_fp) > 0) {
  fp_summary <- state_fp[, .(state, state_label, n_samples)]
  fwrite(fp_summary, file.path(DIRS$reports, "hmm_state_fingerprint_summary.tsv"), sep = "\t")
}

xrf_sig <- xrf_res[q_value < PARAMS$q_threshold]
fwrite(xrf_sig[order(q_value, -abs(beta))], file.path(DIRS$reports, "xrf_significant_associations.tsv"), sep = "\t")

loco_stable <- loco[stability_flag == "stable"]
fwrite(loco_stable[order(q_full, -abs(beta_full))], file.path(DIRS$reports, "loco_stable_associations.tsv"), sep = "\t")

selected_mode <- if (nrow(mode_eval) > 0) mode_eval[1, harmonization_mode] else "NA"
operational_mode <- "corewise_z"

get_origin_n <- function(origin_label) {
  v <- ext_state_origin[state_origin == origin_label, n_samples]
  if (length(v) == 0 || !is.finite(v[1])) return(0L)
  as.integer(v[1])
}

overlap_final_match_pct <- if (nrow(ext_lock_diag) > 0) {
  as.numeric(mean(ext_lock_diag$final_matches_main, na.rm = TRUE) * 100)
} else {
  NA_real_
}

key_metrics <- data.table(
  metric = c(
    "n_samples_main", "n_taxa_main", "n_modules_non_grey",
    "wgcna_input_strategy", "wgcna_profile", "wgcna_training_samples", "hmm_input_balancing_status",
    "preservation_strong", "preservation_moderate", "preservation_weak",
    "hmm_states_main", "r1_r2_state_pct_match", "xrf_significant_n", "loco_stable_n",
    "operational_harmonization_mode", "extended_selected_harmonization_mode", "extended_samples_total", "extended_samples_older_than_150kyr",
    "extended_inherited_main_window_n", "extended_inferred_older_than_150kyr_n", "extended_overlap_final_main_match_pct",
    "extended_scaling_fallback_cores_n",
    "tea_samples_n",
    "tea_module_enrichment_calls_n",
    "tea_taxa_with_oap_n"
  ),
  value = c(
    nsamp, ntaxa, nmods,
    PARAMS$wgcna_input_strategy, PARAMS$wgcna_profile, get_balance_value("total_training_samples"), PARAMS$hmm_input_balancing_status,
    pres_summary$strong, pres_summary$moderate, pres_summary$weak,
    uniqueN(hmm$state),
    concord[metric == "pct_match", value],
    nrow(xrf_sig), nrow(loco_stable),
    operational_mode,
    selected_mode,
    nrow(ext_states),
    nrow(ext_states[age_kyr > PARAMS$main_max_age_kyr]),
    get_origin_n("inherited_main_window"),
    get_origin_n("inferred_extended_only"),
    overlap_final_match_pct,
    uniqueN(ext_scale_diag[fallback_used == TRUE, core]),
    nrow(tea_sample),
    nrow(tea_mod_enrich[analysis_set == "include_grey" & direction %in% c("enriched", "depleted")]),
    tea_taxa_with_oap_n
  )
)
fwrite(key_metrics, file.path(DIRS$reports, "key_metrics.tsv"), sep = "\t")

report_lines <- c(
  "# WGCNA+HMM clean workflow execution report",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Data scope",
  sprintf("- Main-window samples (<= %s kyr): %d", PARAMS$main_max_age_kyr, nsamp),
  sprintf("- Taxa retained after prevalence filter: %d", ntaxa),
  sprintf("- Non-grey modules: %d", nmods),
  sprintf("- WGCNA input strategy: %s", PARAMS$wgcna_input_strategy),
  sprintf("- WGCNA parameter profile: %s", PARAMS$wgcna_profile),
  sprintf("- WGCNA training samples used for module construction: %s", get_balance_value("total_training_samples", "NA")),
  "",
  "## IMPORTANT HMM balancing note",
  "Balanced sampling is currently applied to WGCNA module construction only.",
  "M's HMM stage still consumes all projected main-window module eigengenes.",
  "This HMM input policy must be reviewed with M before final biological claims about HMM states.",
  "",
  "## Module preservation (GeoB25202_R1 -> GeoB25202_R2)",
  sprintf("- Strong: %d", pres_summary$strong),
  sprintf("- Moderate: %d", pres_summary$moderate),
  sprintf("- Weak: %d", pres_summary$weak),
  sprintf("- Median Zsummary: %.3f", pres_summary$median_Zsummary),
  "",
  "## HMM states",
  sprintf("- Operational harmonization mode for states/downstream analyses: %s", operational_mode),
  sprintf("- Number of states in main analysis: %d", uniqueN(hmm$state)),
  sprintf("- R1/R2 interpolated concordance (%% match): %.2f", as.numeric(concord[metric == "pct_match", value])),
  sprintf("- Main-window state fingerprints file present: %s", ifelse(nrow(state_fp) > 0, "yes (state_fingerprints_main.tsv)", "no")),
  "",
  "## XRF and LOCO",
  sprintf("- Significant module-XRF associations (q < %.2f): %d", PARAMS$q_threshold, nrow(xrf_sig)),
  sprintf("- LOCO-stable associations: %d", nrow(loco_stable)),
  "",
  "## Extended-time projection",
  sprintf("- Samples in extended projection: %d", nrow(ext_states)),
  sprintf("- Samples older than %s kyr: %d", PARAMS$main_max_age_kyr, nrow(ext_states[age_kyr > PARAMS$main_max_age_kyr])),
  sprintf("- Locked-overlap policy: states at <= %s kyr are inherited from main analysis; only > %s kyr states are newly inferred", PARAMS$main_max_age_kyr, PARAMS$main_max_age_kyr),
  sprintf("- Inherited main-window states (state_origin=inherited_main_window): %d", get_origin_n("inherited_main_window")),
  sprintf("- Newly inferred older states (state_origin=inferred_extended_only): %d", get_origin_n("inferred_extended_only")),
  sprintf("- Overlap final-state agreement with main (after lock): %s", ifelse(is.finite(overlap_final_match_pct), sprintf("%.2f%%", overlap_final_match_pct), "NA")),
  sprintf("- Selected harmonization mode for extended inference: %s", selected_mode),
  sprintf("- Cores requiring scaling fallback (no core-specific main-window rows): %d", uniqueN(ext_scale_diag[fallback_used == TRUE, core])),
  "- Mode evaluation table: extended_harmonization_mode_evaluation.tsv",
  "- Locked-overlap diagnostics (main vs inferred vs final): extended_locked_overlap_diagnostics.tsv",
  "- State-origin summary (inherited vs inferred): extended_state_origin_summary.tsv",
  "- Scaling parameters used in projection: extended_corewise_z_scaling_parameters_used.tsv",
  "- Scaling fallback diagnostics by core: extended_corewise_z_scaling_fallback_diagnostics.tsv",
  "- Older-period diagnostics by core: extended_older_period_diagnostics.tsv",
  "",
  "## Output index",
  sprintf("- Main results: %s/", DIRS$main),
  sprintf("- Main HMM state fingerprints: %s", file.path(DIRS$main, "state_fingerprints_main.tsv")),
  sprintf("- Main TEA results: %s/", DIRS$main_tea),
  sprintf("- Taxon-level TEA traits: %s", file.path(DIRS$main_tea, "tea_taxon_metabolic_calls.tsv")),
  sprintf("- Extended results: %s/", DIRS$extended),
  sprintf("- Report tables: %s/", DIRS$reports)
)

writeLines(report_lines, con = file.path(DIRS$reports, "workflow_execution_report.md"))

write_run_metadata(
  file.path(DIRS$reports, "06_report_run_metadata.tsv"),
  "06_report.R",
  extra = list(
    q_threshold = PARAMS$q_threshold,
    wgcna_input_strategy = PARAMS$wgcna_input_strategy,
    wgcna_profile = PARAMS$wgcna_profile,
    hmm_input_balancing_status = PARAMS$hmm_input_balancing_status,
    operational_harmonization = operational_mode,
    selected_extended_mode = selected_mode,
    overlap_state_policy = "lock_main_window_states_to_main_analysis",
    overlap_lock_threshold_kyr = PARAMS$main_max_age_kyr,
    overlap_diagnostics_file = "extended_locked_overlap_diagnostics.tsv",
    state_origin_summary_file = "extended_state_origin_summary.tsv"
  )
)
write_session_info(file.path(DIRS$reports, "06_report_sessionInfo.txt"))

log_msg("06_report complete.")
