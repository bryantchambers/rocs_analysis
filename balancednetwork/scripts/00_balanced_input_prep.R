#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(DESeq2)
  library(matrixStats)
})

source(here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

progress_tsv <- Sys.getenv("BALANCED_RUN_PROGRESS_TSV", unset = BALANCED_RUN_PROGRESS_TSV)
if (!nzchar(progress_tsv)) {
  progress_tsv <- file.path(BAL$logs_dir, sprintf("balanced_stage1_prep_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
}

progress_update <- function(stage, status, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      lane = "balancednetwork_stage1_prep",
      stage = stage,
      status = status,
      detail = detail
    ),
    progress_tsv,
    sep = "\t",
    append = file.exists(progress_tsv),
    col.names = !file.exists(progress_tsv)
  )
}

sample_centered_clr <- function(count_mat, pseudocount = 0.5) {
  log_mat <- log(count_mat + pseudocount)
  sweep(log_mat, 2, colMeans(log_mat), "-")
}

force <- identical(Sys.getenv("FORCE", unset = "0"), "1")
metadata_file <- file.path(RESULTS$stage1, "sample_metadata_stage1.tsv")
if (file.exists(metadata_file) && !force) {
  message("[balancednetwork] stage1 prep exists; skipping (use FORCE=1 to rebuild).")
  quit(save = "no")
}

progress_update("stage1_prep", "start", sprintf("variant=%s", BAL$stage1_variant))
log_msg("Loading damage summary...")
tax <- fread(UPSTREAM$tax_damage)
log_msg(sprintf("  %d records, %d columns", nrow(tax), ncol(tax)))

log_msg("Loading metadata...")
meta <- fread(UPSTREAM$metadata)

log_msg("Loading functional classification...")
func <- fread(CLASS$prokaryote_function)

required_cols <- c("label", "subspecies", "is_dmg", "domain", "n_reads", "tax_abund_tad", "tax_abund_read", "reference_length")
missing_cols <- setdiff(required_cols, names(tax))
if (length(missing_cols) > 0L) {
  stop("Damage summary is missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta_s1 <- meta[
  core %in% PARAMS$all_cores &
    y_bp / 1000 <= PARAMS$stage1_max_age_kyr &
    !label %in% PARAMS$excluded_samples
]
log_msg(sprintf("Stage-1 samples (all cores): %d", nrow(meta_s1)))

prok <- tax[
  is_dmg == "Damaged" &
    n_reads >= 100 &
    domain %in% c("d__Archaea", "d__Bacteria", "d__Viruses") &
    label %in% meta_s1$label
]

log_msg(sprintf("  %d rows after damage/read gating", nrow(prok)))

prok_agg <- prok[, .(
  tax_abund_tad = sum(tax_abund_tad),
  tax_abund_read = sum(tax_abund_read),
  reference_length = mean(reference_length)
), by = .(subspecies, label)]

prok_agg[, abund_val := fifelse(tax_abund_tad > 0, tax_abund_tad, tax_abund_read)]
log_msg(sprintf(
  "  %d taxa x %d samples after aggregation",
  uniqueN(prok_agg$subspecies), uniqueN(prok_agg$label)
))

prev_thr <- 10L
keep_taxa <- prok_agg[abund_val > 0, .N, by = subspecies][N >= prev_thr, subspecies]
prok_agg <- prok_agg[subspecies %in% keep_taxa]
log_msg(sprintf("  %d taxa after permissive prevalence filter (>= %d samples)", length(keep_taxa), prev_thr))

wide <- dcast(prok_agg, subspecies ~ label, value.var = "tax_abund_tad", fill = 0L)
taxa_ids <- wide$subspecies
count_mat <- as.matrix(wide[, -1, with = FALSE])
rownames(count_mat) <- taxa_ids

if (any(abs(count_mat - round(count_mat)) > 1e-8, na.rm = TRUE)) {
  stop("`tax_abund_tad` matrix contains non-integer values; DESeq2 input requires count-like integers.")
}
storage.mode(count_mat) <- "integer"

sample_totals_retained <- colSums(count_mat)
zero_samples <- names(sample_totals_retained)[sample_totals_retained <= 0]
if (length(zero_samples) > 0L) {
  log_msg(sprintf("  Removing %d zero-total samples after prevalence filter", length(zero_samples)))
  count_mat <- count_mat[, sample_totals_retained > 0, drop = FALSE]
  sample_totals_retained <- sample_totals_retained[sample_totals_retained > 0]
  meta_s1 <- meta_s1[label %in% colnames(count_mat)]
}

ref_len <- prok_agg[, .(reference_length = mean(reference_length)), by = subspecies]
ref_lengths <- setNames(ref_len[match(taxa_ids, subspecies), reference_length], taxa_ids)

log_msg(sprintf("  Count matrix: %d taxa x %d samples", nrow(count_mat), ncol(count_mat)))
log_msg("CLR transform...")
clr_mat <- sample_centered_clr(count_mat, pseudocount = 0.5)
keep_var <- rowVars(clr_mat) > 0
clr_mat <- clr_mat[keep_var, , drop = FALSE]
log_msg(sprintf("  Removed %d zero-variance taxa", sum(!keep_var)))

vst_mat <- t(clr_mat)
log_msg(sprintf("  Final CLR: %d samples x %d taxa", nrow(vst_mat), ncol(vst_mat)))

log_msg("DESeq2 poscounts normalisation (for EMP, not WGCNA)...")
meta_ord <- meta_s1[match(colnames(count_mat), label)]
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = data.frame(row.names = meta_ord$label, core = factor(meta_ord$core)),
  design = ~1
)
dds <- estimateSizeFactors(dds, type = "poscounts")
len_gm <- exp(mean(log(ref_lengths)))
len_fac <- ref_lengths / len_gm
normalizationFactors(dds) <- outer(len_fac, sizeFactors(dds))

taxa_meta <- func[taxon %in% rownames(clr_mat),
  .(taxon, domain, phylum, functional_group, signal_source, tea_primary, confidence_score, module_completeness_detail)
]

saveRDS(dds, file.path(RESULTS$stage1, "prokaryotes_dds.rds"))
saveRDS(vst_mat, file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
fwrite(taxa_meta, file.path(RESULTS$stage1, "prokaryotes_taxa_metadata.tsv"), sep = "\t")
fwrite(meta_s1, file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"), sep = "\t")

fwrite(
  data.table(
    stage1_variant = BAL$stage1_variant,
    upstream_gate = "is_dmg == 'Damaged' & n_reads >= 100",
    input_measure = "tax_abund_tad",
    prevalence_measure = "abund_val = tax_abund_tad if >0 else tax_abund_read",
    prevalence_rule = sprintf("abund_val > 0 in >= %d samples", prev_thr),
    presence_definition = "abund_val > 0",
    clr_centering = "sample",
    clr_pseudocount = 0.5,
    zero_total_samples_excluded = length(zero_samples),
    matrix_orientation = "samples_by_taxa",
    legacy_matrix_filename = "prokaryotes_vst.rds"
  ),
  file.path(RESULTS$stage1, "wgcna_input_metadata.tsv"),
  sep = "\t"
)

fwrite(
  data.table(
    metric = c(
      "stage1_variant",
      "rows_loaded",
      "rows_after_gate",
      "taxa_after_aggregation",
      "taxa_after_prevalence",
      "samples_retained",
      "zero_total_samples_removed"
    ),
    value = c(
      BAL$stage1_variant,
      nrow(tax),
      nrow(prok),
      uniqueN(prok_agg$subspecies),
      length(keep_taxa),
      ncol(count_mat),
      length(zero_samples)
    )
  ),
  file.path(RESULTS$stage1, "stage1_filter_summary.tsv"),
  sep = "\t"
)

progress_update("stage1_prep", "ok", sprintf("taxa=%d samples=%d", ncol(vst_mat), nrow(vst_mat)))
log_msg(sprintf("Done -> %d samples x %d taxa", nrow(vst_mat), ncol(vst_mat)))
