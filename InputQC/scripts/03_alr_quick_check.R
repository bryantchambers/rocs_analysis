#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(matrixStats)
  library(WGCNA)
  library(ggplot2)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA_character_) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}
n_perm <- as.integer(arg_val("n_perm", "50"))
panel_n <- as.integer(arg_val("panel_n", "20"))

BASE <- here("InputQC")
IN_DIR <- file.path(BASE, "results", "inputs")
OUT <- file.path(BASE, "results", "alr_quick_check")
FIG <- file.path(OUT, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT, "progress.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

eta2_group <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- group[ok]
  if (length(unique(group)) < 2 || var(x) == 0) return(NA_real_)
  grand <- mean(x)
  ss_total <- sum((x - grand)^2)
  if (ss_total == 0) return(NA_real_)
  tab <- data.table(x = x, group = group)
  means <- tab[, .(n = .N, mean_x = mean(x)), by = group]
  sum(means$n * (means$mean_x - grand)^2) / ss_total
}

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

split_by_core <- function(mat, meta) {
  out <- lapply(PARAMS$all_cores, function(core_id) {
    samps <- intersect(meta[core == core_id, label], rownames(mat))
    mat[samps, , drop = FALSE]
  })
  names(out) <- PARAMS$all_cores
  out
}

fit_consensus <- function(expr_by_core, power, deepSplit, mergeCutHeight, minModuleSize) {
  multiExpr <- lapply(PARAMS$stage1_cores, function(core_id) list(data = expr_by_core[[core_id]]))
  names(multiExpr) <- PARAMS$stage1_cores
  blockwiseConsensusModules(
    multiExpr = multiExpr,
    power = power,
    networkType = "signed",
    corType = "pearson",
    maxBlockSize = 5000,
    minModuleSize = minModuleSize,
    deepSplit = deepSplit,
    mergeCutHeight = mergeCutHeight,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )
}

age_aligned_concordance <- function(MEs_all, meta) {
  meta_age <- meta[, .(label, core, age_kyr = y_bp / 1000)]
  dt <- merge(MEs_all, meta_age, by.x = "sample", by.y = "label", all.x = TRUE)
  r1 <- PARAMS$stage1_cores[[length(PARAMS$stage1_cores)]]
  r2 <- PARAMS$validation_core
  me_cols <- setdiff(names(MEs_all), "sample")
  out <- rbindlist(lapply(me_cols, function(me) {
    d1 <- dt[core == r1, .(age_kyr, v = get(me))][order(age_kyr)]
    d2 <- dt[core == r2, .(age_kyr, v = get(me))][order(age_kyr)]
    d1 <- d1[is.finite(age_kyr) & is.finite(v)]
    d2 <- d2[is.finite(age_kyr) & is.finite(v)]
    if (nrow(d1) < 3 || nrow(d2) < 3) {
      return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    }
    lo <- max(min(d1$age_kyr), min(d2$age_kyr))
    hi <- min(max(d1$age_kyr), max(d2$age_kyr))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    }
    x <- seq(lo, hi, length.out = PARAMS$wgcna_stability_age_grid_points)
    y1 <- approx(d1$age_kyr, d1$v, xout = x, rule = 2)$y
    y2 <- approx(d2$age_kyr, d2$v, xout = x, rule = 2)$y
    data.table(
      module = me,
      pearson_r = unname(suppressWarnings(cor(y1, y2, method = "pearson"))),
      spearman_rho = unname(suppressWarnings(cor(y1, y2, method = "spearman"))),
      rmse = sqrt(mean((y1 - y2)^2, na.rm = TRUE))
    )
  }), fill = TRUE)
  setnames(out, old = names(out), new = sub("\\.V1$", "", names(out)))
  out
}

compare_assignments <- function(ref_mods, test_mods) {
  ref <- setNames(ref_mods$module, ref_mods$taxon)
  tst <- setNames(test_mods$module, test_mods$taxon)
  common <- intersect(names(ref), names(tst))
  ref <- ref[common]
  tst <- tst[common]
  ref_bio <- setdiff(sort(unique(ref)), c("grey", "gold"))
  tst_bio <- setdiff(sort(unique(tst)), c("grey", "gold"))

  rows <- rbindlist(lapply(ref_bio, function(m) {
    bg <- names(ref)[ref == m]
    js <- sapply(tst_bio, function(tm) jaccard(bg, names(tst)[tst == tm]))
    if (!length(js)) {
      return(data.table(ref_module = m, best_test_module = NA_character_, best_jaccard = NA_real_, ref_size = length(bg), test_size = NA_integer_))
    }
    best <- names(which.max(js))
    data.table(
      ref_module = m,
      best_test_module = best,
      best_jaccard = unname(max(js, na.rm = TRUE)),
      ref_size = length(bg),
      test_size = sum(tst == best)
    )
  }), fill = TRUE)

  list(
    summary = data.table(
      mean_ref_best_jaccard = mean(rows$best_jaccard, na.rm = TRUE),
      median_ref_best_jaccard = median(rows$best_jaccard, na.rm = TRUE),
      min_ref_best_jaccard = min(rows$best_jaccard, na.rm = TRUE),
      ref_grey_pct = mean(ref == "grey") * 100,
      test_grey_pct = mean(tst == "grey") * 100,
      grey_delta_pct = (mean(tst == "grey") - mean(ref == "grey")) * 100,
      same_grey_fraction = mean((ref == "grey") == (tst == "grey"))
    ),
    details = rows
  )
}

settings <- data.table(
  setting_id = c("exp3", "exp4", "opt5"),
  power = c(12L, 12L, 12L),
  deepSplit = c(3L, 3L, 1L),
  mergeCutHeight = c(0.25, 0.20, 0.20),
  minModuleSize = c(20L, 20L, 20L)
)

count_mat <- readRDS(file.path(IN_DIR, "raw_counts.rds"))
meta <- fread(file.path(IN_DIR, "sample_metadata_stage1.tsv"))
totals <- colSums(count_mat)
log_total <- log10(totals + 1)
log_raw <- log(count_mat + 0.5)
sample_names <- colnames(count_mat)
meta_ord <- meta[match(sample_names, label)]

log_msg("Selecting ALR reference candidates...")
taxa <- rownames(count_mat)
sample_by_taxon <- t(log_raw)
candidate_dt <- data.table(
  taxon = taxa,
  prevalence = rowMeans(count_mat > 0),
  mean_count = rowMeans(count_mat),
  median_count = rowMedians(count_mat),
  var_log = rowVars(log_raw),
  cor_depth = vapply(taxa, function(tx) cor(sample_by_taxon[, tx], log_total, use = "pairwise.complete.obs"), numeric(1)),
  cor_age = vapply(taxa, function(tx) cor(sample_by_taxon[, tx], meta_ord$y_bp, use = "pairwise.complete.obs"), numeric(1)),
  core_eta2 = vapply(taxa, function(tx) eta2_group(sample_by_taxon[, tx], meta_ord$core), numeric(1))
)
candidate_dt[, neutral_score := abs(cor_depth) + abs(cor_age) + fifelse(is.na(core_eta2), 0, core_eta2)]

candidates <- candidate_dt[prevalence >= 0.50 & median_count > 0 & is.finite(neutral_score)]
if (nrow(candidates) < panel_n) {
  candidates <- candidate_dt[prevalence >= 0.25 & mean_count > 0 & is.finite(neutral_score)]
}
setorder(candidates, neutral_score, var_log)
ref_single <- candidates$taxon[1]
ref_panel <- candidates$taxon[seq_len(min(panel_n, nrow(candidates)))]

candidate_dt[, selected := fifelse(taxon == ref_single, "single_reference",
                                   fifelse(taxon %in% ref_panel, "panel_reference", "not_selected"))]
fwrite(candidate_dt[order(neutral_score)], file.path(OUT, "alr_reference_candidates.tsv"), sep = "\t")
fwrite(data.table(reference_type = c("single", "panel"),
                  taxon = c(ref_single, paste(ref_panel, collapse = ",")),
                  n_taxa = c(1L, length(ref_panel))),
       file.path(OUT, "alr_reference_selection.tsv"), sep = "\t")

log_msg(sprintf("Single ALR reference: %s", ref_single))
log_msg(sprintf("Panel ALR reference taxa: %d", length(ref_panel)))

alr_single <- t(sweep(log_raw[setdiff(taxa, ref_single), , drop = FALSE], 2, log_raw[ref_single, ], "-"))
ref_panel_log <- colMeans(log_raw[ref_panel, , drop = FALSE])
alr_panel <- t(sweep(log_raw, 2, ref_panel_log, "-"))

variants <- list(
  alr_single_reference = alr_single,
  alr_panel_reference = alr_panel
)
depth_summary <- rbindlist(lapply(names(variants), function(v) pc_depth_summary(variants[[v]], totals, v)))
fwrite(depth_summary, file.path(OUT, "alr_depth_pca_summary.tsv"), sep = "\t")

summary_rows <- list()
overlap_rows <- list()
pres_rows <- list()
conc_rows <- list()

for (variant_name in names(variants)) {
  mat <- variants[[variant_name]]
  expr_by_core <- split_by_core(mat, meta)
  for (i in seq_len(nrow(settings))) {
    par <- settings[i]
    sid <- par$setting_id
    run_dir <- file.path(OUT, variant_name, sid)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    log_msg(sprintf("Evaluating %s / %s", variant_name, sid))

    net <- fit_consensus(expr_by_core, par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize)
    mods <- data.table(taxon = names(net$colors), module = as.character(net$colors))
    fwrite(mods, file.path(run_dir, "module_assignments.tsv"), sep = "\t")

    ref_mods <- fread(here("networkQC", "results", "full_eval", sid, "module_assignments.tsv"))
    cmp <- compare_assignments(ref_mods, mods)
    details <- cmp$details
    details[, `:=`(variant = variant_name, setting_id = sid)]
    overlap_rows[[length(overlap_rows) + 1]] <- details

    pooled <- do.call(rbind, lapply(PARAMS$stage1_cores, function(core_id) expr_by_core[[core_id]]))
    MEs_train <- orderMEs(moduleEigengenes(pooled, net$colors)$eigengenes)
    MEs_valid <- orderMEs(moduleEigengenes(expr_by_core[[PARAMS$validation_core]], net$colors)$eigengenes)
    MEs_all <- rbind(
      as.data.table(MEs_train, keep.rownames = "sample"),
      as.data.table(MEs_valid, keep.rownames = "sample")
    )
    fwrite(MEs_all, file.path(run_dir, "module_eigengenes.tsv"), sep = "\t")

    var_r1 <- apply(expr_by_core[[PARAMS$stage1_cores[[length(PARAMS$stage1_cores)]]]], 2, var)
    var_r2 <- apply(expr_by_core[[PARAMS$validation_core]], 2, var)
    good <- names(which(var_r1 > 0 & var_r2 > 0))
    mp <- modulePreservation(
      multiData = list(
        R1 = list(data = expr_by_core[[PARAMS$stage1_cores[[length(PARAMS$stage1_cores)]]]][, good]),
        R2 = list(data = expr_by_core[[PARAMS$validation_core]][, good])
      ),
      multiColor = list(R1 = net$colors[good]),
      referenceNetworks = 1,
      testNetworks = 2,
      nPermutations = n_perm,
      randomSeed = PARAMS$seed,
      verbose = 0
    )
    ps <- mp$preservation$Z[[1]][[2]]
    pres <- data.table(
      module = rownames(ps),
      Zsummary = ps$Zsummary.pres,
      Zdensity = ps$Zdensity.pres,
      Zconnectivity = ps$Zconnectivity.pres
    )
    pres[, preserved := fcase(Zsummary > 10, "strong", Zsummary > 2, "moderate", default = "weak")]
    pres[, module_type := fifelse(module %in% c("grey", "gold"), "technical", "biological")]
    pres[, `:=`(variant = variant_name, setting_id = sid)]
    fwrite(pres, file.path(run_dir, "preservation.tsv"), sep = "\t")
    pres_rows[[length(pres_rows) + 1]] <- pres

    conc <- age_aligned_concordance(MEs_all, meta)
    conc[, `:=`(variant = variant_name, setting_id = sid)]
    fwrite(conc, file.path(run_dir, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
    conc_rows[[length(conc_rows) + 1]] <- conc

    bio_pres <- pres[module_type == "biological"]
    bio_conc <- conc[!grepl("grey|gold", module)]
    ds <- depth_summary[variant == variant_name]
    summary_rows[[length(summary_rows) + 1]] <- cbind(
      data.table(
        variant = variant_name,
        setting_id = sid,
        power = par$power,
        deepSplit = par$deepSplit,
        mergeCutHeight = par$mergeCutHeight,
        minModuleSize = par$minModuleSize,
        non_grey_modules = mods[module != "grey", uniqueN(module)],
        grey_pct = mods[module == "grey", .N] / nrow(mods) * 100
      ),
      cmp$summary[, .(
        mean_ref_best_jaccard,
        median_ref_best_jaccard,
        min_ref_best_jaccard,
        grey_delta_pct,
        same_grey_fraction
      )],
      data.table(
        bio_pres_strong = bio_pres[preserved == "strong", .N],
        bio_pres_moderate = bio_pres[preserved == "moderate", .N],
        mean_concordance_pearson = mean(bio_conc$pearson_r, na.rm = TRUE),
        mean_concordance_spearman = mean(bio_conc$spearman_rho, na.rm = TRUE),
        mean_concordance_rmse = mean(bio_conc$rmse, na.rm = TRUE),
        cor_PC1_log_total = ds$cor_PC1_log_total,
        cor_PC2_log_total = ds$cor_PC2_log_total,
        cor_PC3_log_total = ds$cor_PC3_log_total
      )
    )
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
overlap_dt <- rbindlist(overlap_rows, fill = TRUE)
pres_dt <- rbindlist(pres_rows, fill = TRUE)
conc_dt <- rbindlist(conc_rows, fill = TRUE)
summary_dt[, drastic_change := (
  mean_ref_best_jaccard < 0.50 |
    min_ref_best_jaccard < 0.20 |
    abs(grey_delta_pct) > 10 |
    bio_pres_strong < 4 |
    mean_concordance_pearson < 0.60
)]
setorder(summary_dt, setting_id, variant)

fwrite(summary_dt, file.path(OUT, "alr_quick_check_summary.tsv"), sep = "\t")
fwrite(overlap_dt, file.path(OUT, "alr_module_overlap_to_current.tsv"), sep = "\t")
fwrite(pres_dt, file.path(OUT, "alr_preservation_all.tsv"), sep = "\t")
fwrite(conc_dt, file.path(OUT, "alr_eigengene_concordance_all.tsv"), sep = "\t")

heat <- melt(
  summary_dt,
  id.vars = c("variant", "setting_id"),
  measure.vars = c("grey_pct", "mean_ref_best_jaccard", "min_ref_best_jaccard",
                   "bio_pres_strong", "bio_pres_moderate",
                   "mean_concordance_pearson", "mean_concordance_spearman",
                   "mean_concordance_rmse", "cor_PC1_log_total"),
  variable.name = "metric",
  value.name = "value"
)
heat[, scaled := {
  r <- range(value, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) rep(0.5, .N) else (value - r[1]) / (r[2] - r[1])
}, by = metric]
heat[metric %in% c("grey_pct", "mean_concordance_rmse"), scaled := 1 - scaled]
heat[metric == "cor_PC1_log_total", scaled := {
  mx <- max(abs(value), na.rm = TRUE)
  if (!is.finite(mx) || mx == 0) rep(0.5, .N) else 1 - abs(value) / mx
}]
heat[, label := sprintf("%.2f", value)]
heat[, row_label := paste(setting_id, variant, sep = " | ")]
fwrite(heat, file.path(OUT, "alr_quick_check_heatmap_long.tsv"), sep = "\t")

png(file.path(FIG, "alr_quick_check_heatmap.png"), width = 2300, height = 1100, res = 220)
print(
  ggplot(heat, aes(metric, row_label, fill = scaled)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = 2) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f2704", na.value = "grey90", name = "Scaled") +
    labs(
      title = "ALR quick check for WGCNA input sensitivity",
      subtitle = sprintf("Single-taxonomic reference and %d-taxon reference panel; preservation permutations=%d", length(ref_panel), n_perm),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

report <- file.path(BASE, "ALR_QUICK_CHECK_REPORT.md")
sink(report)
cat("# ALR Quick Check Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: ALR-style input sensitivity for `exp3`, `exp4`, and `opt5`\n")
cat("- Preservation permutations: `", n_perm, "`\n\n", sep = "")

cat("## Reference Choice\n\n")
cat("The preprint argues that ALR makes the reference explicit. Here, because we are building a network rather than testing a binary differential-abundance contrast, the reference was chosen as a technical-neutral denominator: high prevalence, nonzero median count, low association with total reads, age, and core.\n\n")
cat("- Single reference taxon: `", ref_single, "`\n", sep = "")
cat("- Reference panel size: `", length(ref_panel), "` taxa\n\n", sep = "")

cat("## Summary\n\n")
cat("|setting|variant|grey_pct|mean_overlap|min_overlap|strong|moderate|pearson|spearman|rmse|PC1_depth|PC2_depth|drastic|\n")
cat("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  cat(sprintf("|%s|%s|%.2f|%.3f|%.3f|%d|%d|%.3f|%.3f|%.3f|%.3f|%.3f|%s|\n",
              r$setting_id, r$variant, r$grey_pct, r$mean_ref_best_jaccard, r$min_ref_best_jaccard,
              r$bio_pres_strong, r$bio_pres_moderate, r$mean_concordance_pearson,
              r$mean_concordance_spearman, r$mean_concordance_rmse,
              r$cor_PC1_log_total, r$cor_PC2_log_total,
              ifelse(isTRUE(r$drastic_change), "yes", "no")))
}
cat("\n")

cat("## Interpretation Guide\n\n")
cat("A good ALR result would reduce depth association across the leading PCs while keeping high overlap with the current best modules, high preservation, and good age-aligned concordance. If ALR strongly changes modules, it is informative but does not automatically solve the preprocessing problem; it means the reference frame controls the network geometry.\n\n")

cat("## Output Map\n\n")
cat("- Reference candidates: `InputQC/results/alr_quick_check/alr_reference_candidates.tsv`\n")
cat("- Summary: `InputQC/results/alr_quick_check/alr_quick_check_summary.tsv`\n")
cat("- Module overlap: `InputQC/results/alr_quick_check/alr_module_overlap_to_current.tsv`\n")
cat("- Preservation: `InputQC/results/alr_quick_check/alr_preservation_all.tsv`\n")
cat("- Concordance: `InputQC/results/alr_quick_check/alr_eigengene_concordance_all.tsv`\n")
cat("- Heatmap: `InputQC/results/alr_quick_check/figures/alr_quick_check_heatmap.png`\n")
sink()

log_msg("ALR quick check complete. Report: ", report)
