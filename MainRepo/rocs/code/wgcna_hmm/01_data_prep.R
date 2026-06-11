#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)
ensure_dirs(c(DIRS$results, DIRS$main, DIRS$extended, DIRS$reports, DIRS$logs))

log_msg("Loading damage table + metadata...")
tax <- if (grepl("\\.gz$", INPUTS$tax_damage)) {
  fread(cmd = sprintf("gzip -dc %s", shQuote(INPUTS$tax_damage)))
} else {
  fread(INPUTS$tax_damage)
}
meta <- fread(INPUTS$metadata)

meta <- meta[
  core %in% PARAMS$all_cores &
    !label %in% PARAMS$excluded_samples
]
meta[, age_kyr := y_bp / 1000]

meta_main <- meta[age_kyr <= PARAMS$main_max_age_kyr]

max_age_map <- data.table(core = names(PARAMS$extended_max_age_kyr), max_age_kyr = as.numeric(PARAMS$extended_max_age_kyr))
meta_ext <- merge(meta, max_age_map, by = "core", all.x = FALSE)
meta_ext <- meta_ext[age_kyr <= max_age_kyr]

log_msg("Filtering to damaged prokaryote/viral taxa for selected samples...")
prok <- tax[
  is_dmg == "Damaged" &
    domain %in% c("d__Archaea", "d__Bacteria", "d__Viruses") &
    label %in% meta_ext$label
]

prok_agg_full <- prok[, .(tax_abund_tad = sum(tax_abund_tad)), by = .(subspecies, label)]

deep_manifest_path <- file.path(DIRS$main, "deep_sediment_excluded_taxa.tsv")
deep_manifest <- data.table(
  taxon = character(),
  exclusion_reason = character(),
  exclusion_detail = character(),
  phylum = character(),
  functional_group = character(),
  display_category = character(),
  ecological_role = character(),
  signal_source = character(),
  confidence_score = numeric(),
  n_samples = integer(),
  total_reads = numeric()
)

n_deep_sediment_candidates <- 0L
n_deep_sediment_excluded_selected <- 0L
n_deep_sediment_excluded_prevalent <- 0L
n_extra_phylum_candidates <- 0L
n_extra_phylum_excluded_selected <- 0L
n_extra_phylum_excluded_prevalent <- 0L

extra_excluded_phyla <- unique(as.character(PARAMS$extra_excluded_phyla))
extra_excluded_phyla <- extra_excluded_phyla[nzchar(extra_excluded_phyla)]

if (isTRUE(PARAMS$exclude_deep_sediment_taxa) || length(extra_excluded_phyla) > 0) {
  log_msg("Sensitivity exclusion filter enabled. Building exclusion set from classification table...")
  class_dt <- fread(INPUTS$prokaryote_function)
  required_cols <- c(
    "taxon", "phylum", "functional_group", "display_category", "ecological_role",
    "signal_source", "confidence_score"
  )
  missing_cols <- setdiff(required_cols, names(class_dt))
  if (length(missing_cols) > 0) {
    stop("Classification table missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  class_dt[, confidence_score := as.numeric(confidence_score)]

  deep_candidates <- data.table(
    taxon = character(), exclusion_reason = character(), exclusion_detail = character(),
    phylum = character(), functional_group = character(), display_category = character(),
    ecological_role = character(), signal_source = character(), confidence_score = numeric()
  )

  if (isTRUE(PARAMS$exclude_deep_sediment_taxa)) {
    deep_candidates <- unique(class_dt[
      display_category == PARAMS$deep_sediment_display_category &
        signal_source %in% PARAMS$deep_sediment_signal_sources &
        confidence_score >= PARAMS$deep_sediment_confidence_min,
      .(
        taxon,
        exclusion_reason = "curated_diagenetic",
        exclusion_detail = sprintf(
          "display_category==%s;signal_source in {%s};confidence_score>=%.2f",
          PARAMS$deep_sediment_display_category,
          paste(PARAMS$deep_sediment_signal_sources, collapse = ","),
          PARAMS$deep_sediment_confidence_min
        ),
        phylum,
        functional_group,
        display_category,
        ecological_role,
        signal_source,
        confidence_score
      )
    ])
    n_deep_sediment_candidates <- uniqueN(deep_candidates$taxon)
  }

  phylum_candidates <- data.table(
    taxon = character(), exclusion_reason = character(), exclusion_detail = character(),
    phylum = character(), functional_group = character(), display_category = character(),
    ecological_role = character(), signal_source = character(), confidence_score = numeric()
  )

  if (length(extra_excluded_phyla) > 0) {
    phylum_candidates <- unique(class_dt[
      phylum %in% extra_excluded_phyla,
      .(
        taxon,
        exclusion_reason = "extra_phylum",
        exclusion_detail = phylum,
        phylum,
        functional_group,
        display_category,
        ecological_role,
        signal_source,
        confidence_score
      )
    ])
    n_extra_phylum_candidates <- uniqueN(phylum_candidates$taxon)
  }

  exclusion_candidates <- unique(rbindlist(list(deep_candidates, phylum_candidates), use.names = TRUE, fill = TRUE))
  taxa_in_run <- unique(prok$subspecies)
  excluded_taxa <- intersect(unique(exclusion_candidates$taxon), taxa_in_run)

  n_deep_sediment_excluded_selected <- length(intersect(unique(deep_candidates$taxon), taxa_in_run))
  n_extra_phylum_excluded_selected <- length(intersect(unique(phylum_candidates$taxon), taxa_in_run))

  if (length(excluded_taxa) > 0) {
    excl_stats <- prok_agg_full[
      subspecies %in% excluded_taxa,
      .(n_samples = sum(tax_abund_tad > 0), total_reads = sum(tax_abund_tad)),
      by = subspecies
    ]

    deep_manifest <- merge(
      exclusion_candidates[taxon %in% excluded_taxa],
      excl_stats[, .(taxon = subspecies, n_samples, total_reads)],
      by = "taxon",
      all.x = TRUE
    )
    deep_manifest[is.na(n_samples), n_samples := 0L]
    deep_manifest[is.na(total_reads), total_reads := 0]

    n_deep_sediment_excluded_prevalent <- prok_agg_full[
      subspecies %in% intersect(unique(deep_candidates$taxon), taxa_in_run) & label %in% meta_main$label & tax_abund_tad > 0,
      .(n_samples = .N),
      by = subspecies
    ][n_samples >= PARAMS$prevalence_min_samples, uniqueN(subspecies)]

    n_extra_phylum_excluded_prevalent <- prok_agg_full[
      subspecies %in% intersect(unique(phylum_candidates$taxon), taxa_in_run) & label %in% meta_main$label & tax_abund_tad > 0,
      .(n_samples = .N),
      by = subspecies
    ][n_samples >= PARAMS$prevalence_min_samples, uniqueN(subspecies)]

    prok <- prok[!subspecies %in% excluded_taxa]
  }

  log_msg(
    "Curated-diagenetic candidates: ", n_deep_sediment_candidates,
    "; excluded in selected samples: ", n_deep_sediment_excluded_selected,
    "; excluded meeting prevalence in main window: ", n_deep_sediment_excluded_prevalent,
    "; extra-phylum candidates: ", n_extra_phylum_candidates,
    "; extra-phylum excluded in selected samples: ", n_extra_phylum_excluded_selected,
    "; extra-phylum excluded meeting prevalence in main window: ", n_extra_phylum_excluded_prevalent
  )
} else {
  log_msg("Sensitivity exclusion filter disabled (default).")
}

setcolorder(
  deep_manifest,
  c("taxon", "exclusion_reason", "exclusion_detail", "phylum", "functional_group", "display_category", "ecological_role", "signal_source", "confidence_score", "n_samples", "total_reads")
)
fwrite(deep_manifest[order(-total_reads, taxon)], deep_manifest_path, sep = "\t")

prok_agg <- prok[, .(tax_abund_tad = sum(tax_abund_tad)), by = .(subspecies, label)]

# prevalence filter on main window to preserve training logic
keep_taxa <- prok_agg[
  label %in% meta_main$label & tax_abund_tad > 0,
  .(n_samples = .N),
  by = subspecies
][n_samples >= PARAMS$prevalence_min_samples, subspecies]

prok_agg <- prok_agg[subspecies %in% keep_taxa]

wide <- dcast(prok_agg, subspecies ~ label, value.var = "tax_abund_tad", fill = 0)
taxa <- wide$subspecies
count_mat <- as.matrix(wide[, -1, with = FALSE])
rownames(count_mat) <- taxa
if (any(abs(count_mat - round(count_mat)) > 1e-8, na.rm = TRUE)) {
  stop("`tax_abund_tad` matrix contains non-integer values; WGCNA-HMM input requires count-like integers.")
}
storage.mode(count_mat) <- "integer"

sample_totals_retained <- colSums(count_mat)
zero_samples <- names(sample_totals_retained)[sample_totals_retained <= 0]
if (length(zero_samples) > 0) {
  log_msg(sprintf("Removing %d zero-total samples after prevalence filter", length(zero_samples)))
  count_mat <- count_mat[, sample_totals_retained > 0, drop = FALSE]
  meta_main <- meta_main[label %in% colnames(count_mat)]
  meta_ext <- meta_ext[label %in% colnames(count_mat)]
}

ext_ids <- intersect(meta_ext$label, colnames(count_mat))

clr_obj <- compute_sample_centered_clr(
  count_mat_taxa_by_samples = count_mat[, ext_ids, drop = FALSE],
  pseudocount = PARAMS$clr_pseudocount
)

clr_mat <- clr_obj$clr

meta_main <- meta_main[label %in% rownames(clr_mat)]
meta_ext <- meta_ext[label %in% rownames(clr_mat)]

write_wgcna_training_design <- function(meta_main, clr_mat) {
  train_meta <- copy(meta_main[core %in% PARAMS$training_cores & label %in% rownames(clr_mat)])
  train_meta[, age_kyr := y_bp / 1000]
  if (nrow(train_meta) == 0) stop("No WGCNA training samples available after filtering.")

  if (identical(PARAMS$wgcna_input_strategy, "original")) {
    manifest <- train_meta[, .(
      label,
      core,
      y_bp,
      age_kyr,
      age_bin = "original_all",
      source = "original_all_training"
    )]
    fwrite(manifest, file.path(DIRS$main, "wgcna_training_samples.tsv"), sep = "\t")
    fwrite(data.table(), file.path(DIRS$main, "balance_bin_counts_long.tsv"), sep = "\t")
    fwrite(data.table(), file.path(DIRS$main, "balance_bin_availability.tsv"), sep = "\t")
    fwrite(data.table(), file.path(DIRS$main, "balance_bin_quotas.tsv"), sep = "\t")
    fwrite(train_meta[, .(label, core, y_bp, age_kyr, selected_for_wgcna = TRUE, reason = "original_all_training")],
           file.path(DIRS$main, "balance_excluded_samples.tsv"), sep = "\t")
    summary_dt <- data.table(
      metric = c("input_strategy", "bin_width_kyr", "retained_bins", "samples_per_core", "total_training_samples", "hmm_input_balancing_status"),
      value = c("original", NA, NA, NA, nrow(manifest), PARAMS$hmm_input_balancing_status)
    )
    fwrite(summary_dt, file.path(DIRS$main, "balance_design_summary.tsv"), sep = "\t")
    return(list(manifest = manifest, summary = summary_dt))
  }

  bin_width_kyr <- PARAMS$wgcna_balance_bin_width_kyr
  if (!is.finite(bin_width_kyr) || bin_width_kyr <= 0) {
    stop("WGCNA_HMM_BALANCE_BIN_WIDTH_KYR must be a positive number.")
  }

  core_range <- train_meta[, .(
    age_min = min(age_kyr, na.rm = TRUE),
    age_max = max(age_kyr, na.rm = TRUE),
    n_available = .N
  ), by = core]
  missing_cores <- setdiff(PARAMS$training_cores, core_range$core)
  if (length(missing_cores) > 0) {
    stop("Missing WGCNA training cores after filtering: ", paste(missing_cores, collapse = ", "))
  }

  shared_age_min <- max(core_range$age_min)
  shared_age_max <- min(core_range$age_max)
  if (!is.finite(shared_age_min) || !is.finite(shared_age_max) || shared_age_max <= shared_age_min) {
    stop("Shared age window could not be computed for balanced WGCNA training.")
  }

  design_pool <- train_meta[age_kyr >= shared_age_min & age_kyr <= shared_age_max]
  design_pool[, age_bin_lo := floor(age_kyr / bin_width_kyr) * bin_width_kyr]
  design_pool[, age_bin_hi := age_bin_lo + bin_width_kyr]
  design_pool[, age_bin := sprintf("%.0f-%.0f", age_bin_lo, age_bin_hi)]

  bin_counts <- design_pool[, .N, by = .(age_bin, age_bin_lo, age_bin_hi, core)]
  setorder(bin_counts, age_bin_lo, core)
  bin_wide <- dcast(bin_counts, age_bin + age_bin_lo + age_bin_hi ~ core, value.var = "N", fill = 0)
  for (core_id in PARAMS$training_cores) {
    if (!core_id %in% names(bin_wide)) bin_wide[, (core_id) := 0L]
  }
  bin_wide[, quota_per_core := do.call(pmin, .SD), .SDcols = PARAMS$training_cores]
  bin_wide[, retained := quota_per_core > 0]
  retained_bins <- bin_wide[retained == TRUE]
  if (nrow(retained_bins) == 0) stop("No balanced age bins with non-zero quota across all training cores.")

  selected <- rbindlist(lapply(seq_len(nrow(retained_bins)), function(i) {
    b <- retained_bins[i]
    rbindlist(lapply(PARAMS$training_cores, function(core_id) {
      pool <- design_pool[core == core_id & age_bin == b$age_bin, .(label, core, y_bp, age_kyr, age_bin)]
      pool[sample.int(.N, size = b$quota_per_core, replace = FALSE)]
    }))
  }))
  selected[, source := "balanced_age_core_quota"]
  setorder(selected, age_kyr, core, label)

  selected_counts <- selected[, .N, by = core][order(core)]
  if (length(unique(selected_counts$N)) != 1L) {
    stop("Balanced WGCNA selection failed core-count equality.")
  }

  all_train <- train_meta[, .(label, core, y_bp, age_kyr)]
  all_train <- merge(all_train, selected[, .(label, selected_for_wgcna = TRUE)], by = "label", all.x = TRUE)
  all_train[is.na(selected_for_wgcna), selected_for_wgcna := FALSE]
  all_train[, reason := fifelse(
    selected_for_wgcna,
    "selected",
    fifelse(age_kyr < shared_age_min | age_kyr > shared_age_max, "outside_shared_age_window", "not_selected_due_to_quota")
  )]

  summary_dt <- data.table(
    metric = c(
      "input_strategy", "bin_width_kyr", "shared_age_min_kyr", "shared_age_max_kyr",
      "retained_bins", "samples_per_core", "total_training_samples", "hmm_input_balancing_status"
    ),
    value = c(
      PARAMS$wgcna_input_strategy,
      bin_width_kyr,
      shared_age_min,
      shared_age_max,
      nrow(retained_bins),
      selected_counts$N[1],
      nrow(selected),
      PARAMS$hmm_input_balancing_status
    )
  )

  fwrite(bin_counts, file.path(DIRS$main, "balance_bin_counts_long.tsv"), sep = "\t")
  fwrite(bin_wide, file.path(DIRS$main, "balance_bin_availability.tsv"), sep = "\t")
  fwrite(retained_bins, file.path(DIRS$main, "balance_bin_quotas.tsv"), sep = "\t")
  fwrite(selected, file.path(DIRS$main, "balanced_baseline_samples.tsv"), sep = "\t")
  fwrite(selected, file.path(DIRS$main, "wgcna_training_samples.tsv"), sep = "\t")
  fwrite(all_train[order(core, age_kyr, label)], file.path(DIRS$main, "balance_excluded_samples.tsv"), sep = "\t")
  fwrite(summary_dt, file.path(DIRS$main, "balance_design_summary.tsv"), sep = "\t")
  fwrite(core_range, file.path(DIRS$main, "balance_core_age_ranges.tsv"), sep = "\t")

  list(manifest = selected, summary = summary_dt)
}

wgcna_design <- write_wgcna_training_design(meta_main, clr_mat)

fwrite(meta_main, file.path(DIRS$main, "sample_metadata_main.tsv"), sep = "\t")
fwrite(meta_ext, file.path(DIRS$extended, "sample_metadata_extended.tsv"), sep = "\t")
fwrite(
  data.table(taxon = colnames(clr_mat)),
  file.path(DIRS$main, "taxa_after_filter.tsv"),
  sep = "\t"
)
saveRDS(clr_mat, file.path(DIRS$main, "clr_matrix_train_centered.rds"))
saveRDS(clr_obj$sample_log_geomeans, file.path(DIRS$main, "clr_sample_log_geomeans.rds"))
fwrite(
  data.table(
    input_measure = "tax_abund_tad",
    prevalence_measure = "tax_abund_tad",
    prevalence_rule = sprintf("tax_abund_tad > 0 in >= %d main-window samples", PARAMS$prevalence_min_samples),
    clr_centering = "sample",
    clr_pseudocount = PARAMS$clr_pseudocount,
    zero_total_samples_excluded = length(zero_samples),
    matrix_orientation = "samples_by_taxa",
    legacy_matrix_filename = "clr_matrix_train_centered.rds"
  ),
  file.path(DIRS$main, "wgcna_input_metadata.tsv"),
  sep = "\t"
)

summary_dt <- data.table(
  metric = c(
    "n_samples_main", "n_samples_extended", "n_taxa_after_prevalence", "zero_total_samples_excluded",
    "wgcna_input_measure", "wgcna_prevalence_measure", "wgcna_clr_centering",
    "wgcna_input_strategy", "wgcna_profile", "wgcna_training_samples", "hmm_input_balancing_status",
    "deep_sediment_filter_enabled", "deep_sediment_candidates_total",
    "deep_sediment_excluded_selected_samples", "deep_sediment_excluded_prevalent_main_window",
    "extra_excluded_phyla_count", "extra_excluded_phyla_values", "extra_phylum_candidates_total",
    "extra_phylum_excluded_selected_samples", "extra_phylum_excluded_prevalent_main_window",
    "main_max_age_kyr", "st13_max_age_kyr_extended", "geob_max_age_kyr_extended"
  ),
  value = c(
    nrow(meta_main), nrow(meta_ext), ncol(clr_mat), length(zero_samples),
    "tax_abund_tad", "tax_abund_tad", "sample",
    PARAMS$wgcna_input_strategy, PARAMS$wgcna_profile, nrow(wgcna_design$manifest), PARAMS$hmm_input_balancing_status,
    as.integer(isTRUE(PARAMS$exclude_deep_sediment_taxa)),
    n_deep_sediment_candidates,
    n_deep_sediment_excluded_selected,
    n_deep_sediment_excluded_prevalent,
    length(extra_excluded_phyla),
    if (length(extra_excluded_phyla) > 0) paste(extra_excluded_phyla, collapse = ";") else "",
    n_extra_phylum_candidates,
    n_extra_phylum_excluded_selected,
    n_extra_phylum_excluded_prevalent,
    PARAMS$main_max_age_kyr,
    PARAMS$extended_max_age_kyr[["ST13"]],
    PARAMS$extended_max_age_kyr[["GeoB25202_R1"]]
  )
)
fwrite(summary_dt, file.path(DIRS$main, "data_prep_summary.tsv"), sep = "\t")

write_run_metadata(
  file.path(DIRS$main, "01_data_prep_run_metadata.tsv"),
  "01_data_prep.R",
  extra = list(
    training_cores = PARAMS$training_cores,
    validation_core = PARAMS$validation_core,
    wgcna_input_strategy = PARAMS$wgcna_input_strategy,
    wgcna_profile = PARAMS$wgcna_profile,
    wgcna_training_samples = nrow(wgcna_design$manifest),
    hmm_input_balancing_status = PARAMS$hmm_input_balancing_status,
    main_max_age_kyr = PARAMS$main_max_age_kyr,
    prevalence_min_samples = PARAMS$prevalence_min_samples,
    clr_pseudocount = PARAMS$clr_pseudocount,
    zero_total_samples_excluded = length(zero_samples),
    wgcna_input_measure = "tax_abund_tad",
    prevalence_measure = "tax_abund_tad",
    clr_centering = "sample",
    exclude_deep_sediment_taxa = PARAMS$exclude_deep_sediment_taxa,
    extra_excluded_phyla = if (length(extra_excluded_phyla) > 0) paste(extra_excluded_phyla, collapse = ";") else "",
    deep_sediment_display_category = PARAMS$deep_sediment_display_category,
    deep_sediment_signal_sources = PARAMS$deep_sediment_signal_sources,
    deep_sediment_confidence_min = PARAMS$deep_sediment_confidence_min,
    deep_sediment_candidates_total = n_deep_sediment_candidates,
    deep_sediment_excluded_selected_samples = n_deep_sediment_excluded_selected,
    deep_sediment_excluded_prevalent_main_window = n_deep_sediment_excluded_prevalent,
    extra_phylum_candidates_total = n_extra_phylum_candidates,
    extra_phylum_excluded_selected_samples = n_extra_phylum_excluded_selected,
    extra_phylum_excluded_prevalent_main_window = n_extra_phylum_excluded_prevalent,
    deep_sediment_manifest = deep_manifest_path
  )
)
write_session_info(file.path(DIRS$main, "01_data_prep_sessionInfo.txt"))

log_msg("01_data_prep complete.")
