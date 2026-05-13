#!/usr/bin/env Rscript
# Copy of scripts/01_data_prep.R adapted for contained NetworkQC input sensitivity.
# It rebuilds the same stage-1 count table, then writes alternative WGCNA inputs.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(DESeq2)
  library(matrixStats)
})

source(here("config.R"))
set.seed(PARAMS$seed)

OUT <- here("networkQC", "input_evaluation", "results", "inputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

pc_depth_summary <- function(mat, totals, variant) {
  common <- intersect(rownames(mat), names(totals))
  pc <- prcomp(mat[common, , drop = FALSE], center = TRUE, scale. = FALSE)$x[, 1:3, drop = FALSE]
  data.table(
    variant = variant,
    n_samples = length(common),
    n_taxa = ncol(mat),
    cor_PC1_log_total = unname(cor(pc[, "PC1"], log10(totals[common] + 1))),
    cor_PC2_log_total = unname(cor(pc[, "PC2"], log10(totals[common] + 1))),
    cor_PC3_log_total = unname(cor(pc[, "PC3"], log10(totals[common] + 1)))
  )
}

log_msg("Loading damage summary, metadata, and functional classification...")
tax <- fread(UPSTREAM$tax_damage)
meta <- fread(UPSTREAM$metadata)
func <- fread(CLASS$prokaryote_function)

meta_s1 <- meta[
  core %in% PARAMS$all_cores &
    y_bp / 1000 <= PARAMS$stage1_max_age_kyr &
    !label %in% PARAMS$excluded_samples
]

prok <- tax[
  is_dmg == "Damaged" &
    domain %in% c("d__Archaea", "d__Bacteria", "d__Viruses") &
    label %in% meta_s1$label
]

prok_agg <- prok[, .(
  n_reads = sum(n_reads),
  reference_length = mean(reference_length)
), by = .(subspecies, label)]

keep_taxa <- prok_agg[n_reads > 0, .(n = .N), by = subspecies][n >= PARAMS$prevalence_min_samples, subspecies]
prok_agg <- prok_agg[subspecies %in% keep_taxa]

wide <- dcast(prok_agg, subspecies ~ label, value.var = "n_reads", fill = 0L)
taxa_ids <- wide$subspecies
count_mat <- as.matrix(wide[, -1, with = FALSE])
rownames(count_mat) <- taxa_ids
storage.mode(count_mat) <- "integer"

ref_len <- prok_agg[, .(reference_length = mean(reference_length)), by = subspecies]
ref_lengths <- setNames(ref_len[match(taxa_ids, subspecies), reference_length], taxa_ids)

log_msg(sprintf("Stage-1 count matrix: %d taxa x %d samples", nrow(count_mat), ncol(count_mat)))

meta_ord <- meta_s1[match(colnames(count_mat), label)]
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = data.frame(row.names = meta_ord$label, core = factor(meta_ord$core)),
  design = ~ 1
)
dds <- estimateSizeFactors(dds, type = "poscounts")
len_gm <- exp(mean(log(ref_lengths)))
len_fac <- ref_lengths / len_gm
normalizationFactors(dds) <- outer(len_fac, sizeFactors(dds))

totals <- colSums(count_mat)
log_total <- log10(totals + 1)

# Current main-pipeline input: taxon-centered log raw counts.
log_raw <- log(count_mat + 0.5)
current_taxon_centered <- t(log_raw - rowMeans(log_raw))

# Standard sample-wise CLR on raw filtered counts.
sample_clr_raw <- t(sweep(log_raw, 2, colMeans(log_raw), "-"))

# DESeq2 size-factor + reference-length normalized log abundance, taxon-centered.
norm_counts <- counts(dds, normalized = TRUE)
log_norm <- log(norm_counts + 0.5)
deseq_length_log <- t(log_norm - rowMeans(log_norm))

# Remove linear log-depth signal per taxon from the current log matrix.
sample_by_taxon <- t(log_raw)
x <- cbind(intercept = 1, log_total = as.numeric(log_total[colnames(count_mat)]))
coef <- solve(crossprod(x), crossprod(x, sample_by_taxon))
log_depth_residualized <- sample_by_taxon - x %*% coef
rownames(log_depth_residualized) <- colnames(count_mat)
colnames(log_depth_residualized) <- rownames(count_mat)

variants <- list(
  current_taxon_centered_log = current_taxon_centered,
  sample_clr_raw = sample_clr_raw,
  deseq_length_log = deseq_length_log,
  log_depth_residualized = log_depth_residualized
)

variant_summary <- rbindlist(lapply(names(variants), function(v) {
  mat <- variants[[v]]
  keep_var <- colVars(mat) > 0
  mat <- mat[, keep_var, drop = FALSE]
  variants[[v]] <<- mat
  saveRDS(mat, file.path(OUT, paste0(v, ".rds")))
  pc_depth_summary(mat, totals, v)
}))

taxa_meta <- func[taxon %in% taxa_ids,
                  .(taxon, domain, phylum, functional_group,
                    signal_source, tea_primary, confidence_score,
                    module_completeness_detail)]

saveRDS(count_mat, file.path(OUT, "raw_counts.rds"))
saveRDS(dds, file.path(OUT, "prokaryotes_dds_input_eval.rds"))
fwrite(meta_s1, file.path(OUT, "sample_metadata_stage1.tsv"), sep = "\t")
fwrite(taxa_meta, file.path(OUT, "prokaryotes_taxa_metadata.tsv"), sep = "\t")
fwrite(data.table(sample = names(totals), total_reads = as.numeric(totals), log_total_reads = log_total),
       file.path(OUT, "sample_depth_summary.tsv"), sep = "\t")
fwrite(variant_summary, file.path(OUT, "input_variant_depth_pca_summary.tsv"), sep = "\t")

log_msg("Input variants written to ", OUT)
