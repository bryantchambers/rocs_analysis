library(dplyr)
library(readr)
library(tibble)

ROCS_ROOT <- Sys.getenv("ROCS_ROOT", "/maps/projects/caeg/people/ngm902/apps/repos/rocs")
ROCS_ROOT <- normalizePath(ROCS_ROOT, mustWork = FALSE)

input_file <- file.path(ROCS_ROOT, "results/microbial/taxonomy/dmg-summary-ssp_clean.tsv.gz")
out_dir <- "/projects/caeg/people/ngm902/apps/repos/rocs/results/microbial/damage/damage-classification-depositional"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")

log_msg <- function(...) {
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = log_file, append = TRUE)
}

clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

relu <- function(x) {
  pmax(0, x)
}

params <- list(
  model_version = "localfixed-depositional-ecofree-v1",
  local_ab_threshold = 0.1,
  local_zfit_threshold = 2,
  sample_norm_eps_ab = 0.02,
  sample_norm_eps_zfit = 0.25,
  sample_norm_eps_rho = 0.10,
  eb_prior_strength = 20,
  score_intercept = -1.20,
  score_w_ref_logit = 1.30,
  score_w_ab_rel = 0.95,
  score_w_zfit_rel = 0.95,
  score_w_rho_rel = 0.40,
  score_w_lowdamage_interaction = 0.55,
  depositional_non_damaged_threshold = 0.33,
  depositional_damaged_threshold_base = 0.75,
  depositional_damaged_threshold_relax_max = 0.05,
  random_seed = 20260324
)

set.seed(params$random_seed)

required_cols <- c("taxid", "label", "A_b", "Zfit", "rho_c")

log_msg("Reading input: ", input_file)
dat <- read_tsv(input_file, show_col_types = FALSE)
dat_input <- dat %>%
  mutate(.row_id = row_number())

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat %>%
  mutate(
    .row_id = row_number(),
    taxid = as.character(taxid),
    A_b = as.numeric(A_b),
    Zfit = as.numeric(Zfit),
    rho_c_adj = ifelse(is.na(rho_c), 0, as.numeric(rho_c)),
    Zfit_adj = ifelse(is.na(Zfit), 0, as.numeric(Zfit)),
    local_fixed_success = as.integer(A_b > params$local_ab_threshold & Zfit_adj >= params$local_zfit_threshold),
    local_fixed_label = ifelse(local_fixed_success == 1, "Damaged", "Non-damaged")
  )

global_local_rate <- mean(dat$local_fixed_success, na.rm = TRUE)

sample_baseline <- dat %>%
  group_by(label) %>%
  summarise(
    sample_A_b_median = median(A_b, na.rm = TRUE),
    sample_A_b_iqr = IQR(A_b, na.rm = TRUE),
    sample_Zfit_median = median(Zfit_adj, na.rm = TRUE),
    sample_Zfit_iqr = IQR(Zfit_adj, na.rm = TRUE),
    sample_rho_median = median(rho_c_adj, na.rm = TRUE),
    sample_rho_iqr = IQR(rho_c_adj, na.rm = TRUE),
    sample_local_fixed_rate = mean(local_fixed_success, na.rm = TRUE),
    n_rows_sample = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sample_A_b_iqr = ifelse(is.na(sample_A_b_iqr), 0, sample_A_b_iqr),
    sample_Zfit_iqr = ifelse(is.na(sample_Zfit_iqr), 0, sample_Zfit_iqr),
    sample_rho_iqr = ifelse(is.na(sample_rho_iqr), 0, sample_rho_iqr),
    low_damage_index = ifelse(
      global_local_rate > 0,
      clamp((global_local_rate - sample_local_fixed_rate) / (global_local_rate + 1e-6), 0, 1),
      1
    )
  )

p0 <- global_local_rate

reference_recurrence <- dat %>%
  group_by(taxid) %>%
  summarise(
    n_ref = n(),
    n_samples_ref = n_distinct(label),
    n_cores_ref = if ("core" %in% names(dat)) n_distinct(core) else NA_integer_,
    local_fixed_successes = sum(local_fixed_success, na.rm = TRUE),
    local_fixed_failures = n_ref - local_fixed_successes,
    local_fixed_rate_empirical = mean(local_fixed_success, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    recurrence_prior_mean = p0,
    recurrence_posterior = (local_fixed_successes + params$eb_prior_strength * recurrence_prior_mean) /
      (n_ref + params$eb_prior_strength),
    recurrence_logit = qlogis(clamp(recurrence_posterior, 1e-6, 1 - 1e-6)),
    recurrence_support_weight = n_ref / (n_ref + params$eb_prior_strength)
  )

dat_scored <- dat %>%
  left_join(sample_baseline, by = "label") %>%
  left_join(reference_recurrence %>% select(taxid, recurrence_posterior, recurrence_logit, recurrence_support_weight, n_ref, local_fixed_successes, local_fixed_rate_empirical), by = "taxid") %>%
  mutate(
    A_b_rel = (A_b - sample_A_b_median) / (sample_A_b_iqr + params$sample_norm_eps_ab),
    Zfit_rel = (Zfit_adj - sample_Zfit_median) / (sample_Zfit_iqr + params$sample_norm_eps_zfit),
    rho_rel = (rho_c_adj - sample_rho_median) / (sample_rho_iqr + params$sample_norm_eps_rho),
    A_b_rel = ifelse(is.finite(A_b_rel), A_b_rel, 0),
    Zfit_rel = ifelse(is.finite(Zfit_rel), Zfit_rel, 0),
    rho_rel = ifelse(is.finite(rho_rel), rho_rel, 0),
    recurrence_posterior = ifelse(is.na(recurrence_posterior), p0, recurrence_posterior),
    recurrence_logit = ifelse(is.na(recurrence_logit), qlogis(clamp(p0, 1e-6, 1 - 1e-6)), recurrence_logit),
    recurrence_support_weight = ifelse(is.na(recurrence_support_weight), 0, recurrence_support_weight),
    n_ref = ifelse(is.na(n_ref), 0, n_ref),
    local_fixed_successes = ifelse(is.na(local_fixed_successes), 0, local_fixed_successes),
    local_fixed_rate_empirical = ifelse(is.na(local_fixed_rate_empirical), 0, local_fixed_rate_empirical),
    depositional_damage_eligible = local_fixed_rate_empirical >= 0.75 & n_ref >= 8,
    relative_row_evidence =
      params$score_w_ab_rel * relu(A_b_rel) +
      params$score_w_zfit_rel * relu(Zfit_rel) +
      params$score_w_rho_rel * relu(rho_rel),
    depositional_score_raw =
      params$score_intercept +
      params$score_w_ref_logit * recurrence_logit +
      relative_row_evidence +
      params$score_w_lowdamage_interaction * low_damage_index * relative_row_evidence,
    depositional_prob_ecofree = inv_logit(depositional_score_raw),
    depositional_damaged_threshold = pmax(
      params$depositional_non_damaged_threshold + 0.05,
      params$depositional_damaged_threshold_base -
        params$depositional_damaged_threshold_relax_max * low_damage_index
    ),
    depositional_label_ecofree = case_when(
      depositional_prob_ecofree >= depositional_damaged_threshold & depositional_damage_eligible ~ "Damaged",
      depositional_prob_ecofree <= params$depositional_non_damaged_threshold ~ "Non-damaged",
      TRUE ~ "Uncertain"
    ),
    downstream_binary_label = ifelse(
      local_fixed_label == "Damaged" | depositional_label_ecofree == "Damaged",
      "Damaged",
      "Non-damaged"
    ),
    rescued_by_deposition = local_fixed_label == "Non-damaged" & depositional_label_ecofree == "Damaged"
  )

taxonomy_cols <- c(
  "taxid", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "name"
)
taxonomy_cols <- taxonomy_cols[taxonomy_cols %in% names(dat_scored)]

reference_summary <- dat_scored %>%
  group_by(taxid) %>%
  summarise(
    n_observations = n(),
    n_samples = n_distinct(label),
    n_cores = if ("core" %in% names(dat_scored)) n_distinct(core) else NA_integer_,
    local_fixed_successes = sum(local_fixed_success, na.rm = TRUE),
    local_fixed_rate = mean(local_fixed_success, na.rm = TRUE),
    recurrence_posterior = first(recurrence_posterior),
    recurrence_support_weight = first(recurrence_support_weight),
    recurrence_logit = first(recurrence_logit),
    mean_A_b_rel = mean(A_b_rel, na.rm = TRUE),
    mean_Zfit_rel = mean(Zfit_rel, na.rm = TRUE),
    mean_rho_rel = mean(rho_rel, na.rm = TRUE),
    mean_depositional_prob_ecofree = mean(depositional_prob_ecofree, na.rm = TRUE),
    median_depositional_prob_ecofree = median(depositional_prob_ecofree, na.rm = TRUE),
    frac_depositional_damaged = mean(depositional_label_ecofree == "Damaged", na.rm = TRUE),
    frac_depositional_uncertain = mean(depositional_label_ecofree == "Uncertain", na.rm = TRUE),
    frac_depositional_non_damaged = mean(depositional_label_ecofree == "Non-damaged", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(dat_scored %>% distinct(across(all_of(taxonomy_cols))), by = "taxid") %>%
  select(all_of(taxonomy_cols), everything()) %>%
  arrange(desc(recurrence_posterior), desc(n_observations))

label_crosstab <- dat_scored %>%
  count(local_fixed_label, depositional_label_ecofree, downstream_binary_label, name = "n") %>%
  mutate(frac_total = n / sum(n)) %>%
  arrange(desc(n), local_fixed_label, depositional_label_ecofree)

topline_summary <- tibble(
  n_rows = nrow(dat_scored),
  n_samples = n_distinct(dat_scored$label),
  n_taxa = n_distinct(dat_scored$taxid),
  global_local_fixed_rate = global_local_rate,
  n_local_fixed_damaged = sum(dat_scored$local_fixed_label == "Damaged", na.rm = TRUE),
  n_depositional_damaged = sum(dat_scored$depositional_label_ecofree == "Damaged", na.rm = TRUE),
  n_depositional_uncertain = sum(dat_scored$depositional_label_ecofree == "Uncertain", na.rm = TRUE),
  n_depositional_non_damaged = sum(dat_scored$depositional_label_ecofree == "Non-damaged", na.rm = TRUE),
  n_rescued_by_deposition = sum(dat_scored$rescued_by_deposition, na.rm = TRUE),
  frac_rescued_of_local_non_damaged = mean(dat_scored$rescued_by_deposition[dat_scored$local_fixed_label == "Non-damaged"], na.rm = TRUE),
  n_downstream_damaged = sum(dat_scored$downstream_binary_label == "Damaged", na.rm = TRUE),
  n_downstream_non_damaged = sum(dat_scored$downstream_binary_label == "Non-damaged", na.rm = TRUE)
)

metadata <- tibble(
  model_version = params$model_version,
  run_timestamp = as.character(Sys.time()),
  input_file = input_file,
  input_md5 = unname(tools::md5sum(input_file)),
  local_rule = sprintf("Damaged if A_b > %.3f and Zfit >= %.3f; else Non-damaged", params$local_ab_threshold, params$local_zfit_threshold),
  depositional_rule = "Ecology-free fixed-score model using recurrence posterior + sample-normalized molecular evidence",
  recurrence_definition = sprintf("Posterior = (successes + %.1f * p0) / (n_ref + %.1f), where successes are local-fixed Damaged calls only", params$eb_prior_strength, params$eb_prior_strength),
  depositional_probability_formula = "inv_logit(intercept + w_ref*recurrence_logit + row_evidence + w_low*low_damage_index*row_evidence)",
  row_evidence_formula = "w_ab*relu(A_b_rel) + w_zfit*relu(Zfit_rel) + w_rho*relu(rho_rel)",
  score_intercept = params$score_intercept,
  score_w_ref_logit = params$score_w_ref_logit,
  score_w_ab_rel = params$score_w_ab_rel,
  score_w_zfit_rel = params$score_w_zfit_rel,
  score_w_rho_rel = params$score_w_rho_rel,
  score_w_lowdamage_interaction = params$score_w_lowdamage_interaction,
  depositional_non_damaged_threshold = params$depositional_non_damaged_threshold,
  depositional_damaged_threshold_base = params$depositional_damaged_threshold_base,
  depositional_damaged_threshold_relax_max = params$depositional_damaged_threshold_relax_max,
  depositional_damage_eligibility_rule = "Reference must have local_fixed_successes / n_observations >= 0.75 and n_observations >= 8",
  random_seed = params$random_seed,
  n_rows = nrow(dat_scored),
  n_samples = n_distinct(dat_scored$label),
  n_taxa = n_distinct(dat_scored$taxid)
)

row_output <- dat_scored %>%
  transmute(
    .row_id,
    is_dmg_local = local_fixed_label,
    sample_A_b_median,
    sample_A_b_iqr,
    depositional_prob = depositional_prob_ecofree,
    depositional_threshold = depositional_damaged_threshold,
    depositional_label = depositional_label_ecofree,
    is_dmg = ifelse(local_fixed_label == "Damaged" | depositional_label_ecofree == "Damaged", "Damaged", "Non-damaged")
  )

main_output <- dat_input %>%
  left_join(row_output, by = ".row_id") %>%
  select(-.row_id)

details_output <- dat_scored %>%
  select(.row_id, all_of(setdiff(names(dat_scored), names(dat_input)))) %>%
  arrange(.row_id)

row_file <- file.path(out_dir, "dmg-summary-ssp-damage-classification-depositional.tsv.gz")
details_file <- file.path(out_dir, "dmg-summary-ssp-damage-classification-depositional-details.tsv.gz")
sample_baseline_file <- file.path(out_dir, "sample-baselines.tsv")
ref_file <- file.path(out_dir, "reference-summary.tsv")
crosstab_file <- file.path(out_dir, "label-crosstab.tsv")
summary_file <- file.path(out_dir, "summary.tsv")
metadata_file <- file.path(out_dir, "metadata.tsv")
session_file <- file.path(out_dir, "sessioninfo.txt")

log_msg("Writing outputs to: ", out_dir)
write_tsv(main_output, row_file)
write_tsv(details_output, details_file)
write_tsv(sample_baseline, sample_baseline_file)
write_tsv(reference_summary, ref_file)
write_tsv(label_crosstab, crosstab_file)
write_tsv(topline_summary, summary_file)
write_tsv(metadata, metadata_file)
writeLines(capture.output(sessionInfo()), session_file)

checksum_targets <- c(
  row_file,
  details_file,
  sample_baseline_file,
  ref_file,
  crosstab_file,
  summary_file,
  metadata_file,
  session_file
)

checksums <- tibble(
  file = basename(checksum_targets),
  md5 = unname(tools::md5sum(checksum_targets))
)

checksums_file <- file.path(out_dir, "output-checksums.tsv")
write_tsv(checksums, checksums_file)

log_msg("Completed successfully")

