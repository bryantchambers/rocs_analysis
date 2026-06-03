#!/usr/bin/env Rscript
# 01_data_prep.R — Filter taxa, normalise, CLR-transform, build trait table
#
# Input:   dmg-summary-ssp_selected.tsv.gz  (pre-filtered to damaged reads,
#          is_dmg == "Damaged" already set by upstream CCC criterion)
# Method:  estimateSizeFactors(type="poscounts") + reference-length offset → CLR
# Outputs: results/stage1/
#   prokaryotes_dds.rds           DESeq2 object (used by 04_emp.R for norm counts)
#   prokaryotes_vst.rds           CLR matrix (samples × taxa); named "vst" for legacy
#   prokaryotes_taxa_metadata.tsv per-taxon functional annotations
#   sample_metadata_stage1.tsv    filtered sample table (all 4 cores ≤ 150 kyr)

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(DESeq2)
  library(matrixStats)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── 1. Load ───────────────────────────────────────────────────────────────────

log_msg("Loading damage summary...")
tax <- fread(UPSTREAM$tax_damage)
log_msg(sprintf("  %d records, %d columns", nrow(tax), ncol(tax)))

log_msg("Loading metadata...")
meta <- fread(UPSTREAM$metadata)

log_msg("Loading functional classification...")
func <- fread(CLASS$prokaryote_function)

# ── 2. Stage-1 sample filter ──────────────────────────────────────────────────

meta_s1 <- meta[
  core %in% PARAMS$all_cores &
  y_bp / 1000 <= PARAMS$stage1_max_age_kyr &
  !label %in% PARAMS$excluded_samples
]
log_msg(sprintf("Stage-1 samples (all cores): %d", nrow(meta_s1)))

# ── 3. Filter records ─────────────────────────────────────────────────────────

prok <- tax[
  is_dmg == "Damaged" &
  domain %in% c("d__Archaea", "d__Bacteria", "d__Viruses") &
  label %in% meta_s1$label
]

# Aggregate by subspecies × sample
prok_agg <- prok[, .(
  n_reads          = sum(n_reads),
  reference_length = mean(reference_length)
), by = .(subspecies, label)]

log_msg(sprintf("  %d taxa × %d samples (after damage + domain filter)",
                uniqueN(prok_agg$subspecies), uniqueN(prok_agg$label)))

# ── 4. Prevalence filter (≥10 samples, absolute count — matches original) ────

prev_thr  <- PARAMS$prevalence_min_samples
keep_taxa <- prok_agg[n_reads > 0, .(n = .N), by = subspecies][n >= prev_thr, subspecies]
prok_agg  <- prok_agg[subspecies %in% keep_taxa]

log_msg(sprintf("  %d taxa after prevalence filter (>= %d samples)", length(keep_taxa), prev_thr))

# ── 5. Count matrix (taxa × samples) ─────────────────────────────────────────

wide      <- dcast(prok_agg, subspecies ~ label, value.var = "n_reads", fill = 0L)
taxa_ids  <- wide$subspecies
count_mat <- as.matrix(wide[, -1, with = FALSE])
rownames(count_mat) <- taxa_ids
storage.mode(count_mat) <- "integer"

# Reference lengths
ref_len <- prok_agg[, .(reference_length = mean(reference_length)), by = subspecies]
ref_lengths <- setNames(ref_len[match(taxa_ids, subspecies), reference_length], taxa_ids)

log_msg(sprintf("  Count matrix: %d taxa × %d samples", nrow(count_mat), ncol(count_mat)))

# ── 6. CLR transform directly on raw counts (matches original) ───────────────
# The original applies CLR to raw filtered counts, not to DESeq2-normalised
# counts. CLR itself handles compositionality, so library-size normalisation
# before CLR is redundant and would alter the result.

log_msg("CLR transform...")

# Pseudocount 0.5 before log: prevents log(0) with minimal distortion of
# abundant taxa; matches original (clr_pseudocount = 0.5).
# Subtracting row means centres each taxon across samples (geometric mean
# ratio transform), making the matrix compositionally coherent for WGCNA.
clr_mat <- log(count_mat + 0.5)
clr_mat <- clr_mat - rowMeans(clr_mat)   # taxa × samples

# Remove zero-variance taxa
keep_var <- rowVars(clr_mat) > 0
clr_mat  <- clr_mat[keep_var, ]
log_msg(sprintf("  Removed %d zero-variance taxa", sum(!keep_var)))

# Transpose to samples × taxa (WGCNA convention)
vst_mat <- t(clr_mat)
log_msg(sprintf("  Final CLR: %d samples × %d taxa", nrow(vst_mat), ncol(vst_mat)))

# ── 7. DESeq2 poscounts + reference-length normalisation (for EMP in 04_emp.R) ──
# Built separately from WGCNA input; the original also builds DESeq2 here for
# downstream projection/EMP use, not for the WGCNA matrix.

log_msg("DESeq2 poscounts normalisation (for EMP, not WGCNA)...")

meta_ord <- meta_s1[match(colnames(count_mat), label)]

# design = ~ 1: intercept-only unconditional size factors.
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = data.frame(row.names = meta_ord$label, core = factor(meta_ord$core)),
  design    = ~ 1
)

# poscounts handles zero-inflated ancient DNA
dds <- estimateSizeFactors(dds, type = "poscounts")

# Combined normalisation factor: size factor × reference-length factor.
# Reference lengths vary ~10× across taxa; dividing by the geometric mean length
# puts taxa on a common per-bp scale before multiplying by the poscounts size
# factor. outer() produces a taxa × samples matrix for DESeq2's slot convention.
len_gm  <- exp(mean(log(ref_lengths)))
len_fac <- ref_lengths / len_gm
nf      <- outer(len_fac, sizeFactors(dds))   # taxa × samples
normalizationFactors(dds) <- nf

# ── 8. Per-taxon metadata ────────────────────────────────────────────────────

taxa_meta <- func[taxon %in% rownames(clr_mat),
                   .(taxon, domain, phylum, functional_group,
                     signal_source, tea_primary, confidence_score,
                     module_completeness_detail)]

# ── 9. Save ───────────────────────────────────────────────────────────────────

log_msg("Saving...")
saveRDS(dds,     file.path(RESULTS$stage1, "prokaryotes_dds.rds"))
saveRDS(vst_mat, file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
fwrite(taxa_meta, file.path(RESULTS$stage1, "prokaryotes_taxa_metadata.tsv"), sep = "\t")
fwrite(meta_s1,   file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"),    sep = "\t")

log_msg(sprintf("Done → %d samples × %d taxa", nrow(vst_mat), ncol(vst_mat)))
