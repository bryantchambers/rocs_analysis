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

prok_agg_full <- prok[, .(n_reads = sum(n_reads)), by = .(subspecies, label)]

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
      .(n_samples = sum(n_reads > 0), total_reads = sum(n_reads)),
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
      subspecies %in% intersect(unique(deep_candidates$taxon), taxa_in_run) & label %in% meta_main$label & n_reads > 0,
      .(n_samples = .N),
      by = subspecies
    ][n_samples >= PARAMS$prevalence_min_samples, uniqueN(subspecies)]

    n_extra_phylum_excluded_prevalent <- prok_agg_full[
      subspecies %in% intersect(unique(phylum_candidates$taxon), taxa_in_run) & label %in% meta_main$label & n_reads > 0,
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

prok_agg <- prok[, .(n_reads = sum(n_reads)), by = .(subspecies, label)]

# prevalence filter on main window to preserve training logic
keep_taxa <- prok_agg[
  label %in% meta_main$label & n_reads > 0,
  .(n_samples = .N),
  by = subspecies
][n_samples >= PARAMS$prevalence_min_samples, subspecies]

prok_agg <- prok_agg[subspecies %in% keep_taxa]

wide <- dcast(prok_agg, subspecies ~ label, value.var = "n_reads", fill = 0)
taxa <- wide$subspecies
count_mat <- as.matrix(wide[, -1, with = FALSE])
rownames(count_mat) <- taxa
storage.mode(count_mat) <- "integer"

main_ids <- intersect(meta_main$label, colnames(count_mat))
ext_ids <- intersect(meta_ext$label, colnames(count_mat))

clr_obj <- compute_clr_train_centered(
  count_mat_taxa_by_samples = count_mat[, ext_ids, drop = FALSE],
  train_sample_ids = main_ids,
  pseudocount = PARAMS$clr_pseudocount
)

clr_mat <- clr_obj$clr

meta_main <- meta_main[label %in% rownames(clr_mat)]
meta_ext <- meta_ext[label %in% rownames(clr_mat)]

fwrite(meta_main, file.path(DIRS$main, "sample_metadata_main.tsv"), sep = "\t")
fwrite(meta_ext, file.path(DIRS$extended, "sample_metadata_extended.tsv"), sep = "\t")
fwrite(
  data.table(taxon = colnames(clr_mat)),
  file.path(DIRS$main, "taxa_after_filter.tsv"),
  sep = "\t"
)
saveRDS(clr_mat, file.path(DIRS$main, "clr_matrix_train_centered.rds"))
saveRDS(clr_obj$train_row_means, file.path(DIRS$main, "clr_training_row_means.rds"))

summary_dt <- data.table(
  metric = c(
    "n_samples_main", "n_samples_extended", "n_taxa_after_prevalence",
    "deep_sediment_filter_enabled", "deep_sediment_candidates_total",
    "deep_sediment_excluded_selected_samples", "deep_sediment_excluded_prevalent_main_window",
    "extra_excluded_phyla_count", "extra_excluded_phyla_values", "extra_phylum_candidates_total",
    "extra_phylum_excluded_selected_samples", "extra_phylum_excluded_prevalent_main_window",
    "main_max_age_kyr", "st13_max_age_kyr_extended", "geob_max_age_kyr_extended"
  ),
  value = c(
    nrow(meta_main), nrow(meta_ext), ncol(clr_mat),
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
    main_max_age_kyr = PARAMS$main_max_age_kyr,
    prevalence_min_samples = PARAMS$prevalence_min_samples,
    clr_pseudocount = PARAMS$clr_pseudocount,
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
