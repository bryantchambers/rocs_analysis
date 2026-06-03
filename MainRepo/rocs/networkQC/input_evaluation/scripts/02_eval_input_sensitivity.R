#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
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

n_perm <- as.integer(arg_val("n_perm", "100"))

BASE <- here("networkQC", "input_evaluation")
IN_DIR <- file.path(BASE, "results", "inputs")
OUT <- file.path(BASE, "results", "sensitivity")
FIG <- file.path(BASE, "results", "figures")
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

  ref_rows <- rbindlist(lapply(ref_bio, function(m) {
    bg <- names(ref)[ref == m]
    js <- sapply(tst_bio, function(tm) jaccard(bg, names(tst)[tst == tm]))
    if (!length(js)) return(data.table(ref_module = m, best_test_module = NA_character_, best_jaccard = NA_real_, ref_size = length(bg), test_size = NA_integer_))
    best <- names(which.max(js))
    data.table(ref_module = m, best_test_module = best, best_jaccard = unname(max(js, na.rm = TRUE)),
               ref_size = length(bg), test_size = sum(tst == best))
  }), fill = TRUE)

  data.table(
    mean_ref_best_jaccard = mean(ref_rows$best_jaccard, na.rm = TRUE),
    median_ref_best_jaccard = median(ref_rows$best_jaccard, na.rm = TRUE),
    min_ref_best_jaccard = min(ref_rows$best_jaccard, na.rm = TRUE),
    ref_grey_pct = mean(ref == "grey") * 100,
    test_grey_pct = mean(tst == "grey") * 100,
    grey_delta_pct = (mean(tst == "grey") - mean(ref == "grey")) * 100,
    same_grey_fraction = mean((ref == "grey") == (tst == "grey")),
    ref_module_rows = list(ref_rows)
  )
}

settings <- data.table(
  setting_id = c("exp3", "exp4", "opt5"),
  power = c(12L, 12L, 12L),
  deepSplit = c(3L, 3L, 1L),
  mergeCutHeight = c(0.25, 0.20, 0.20),
  minModuleSize = c(20L, 20L, 20L)
)

variants <- c("current_taxon_centered_log", "sample_clr_raw", "deseq_length_log", "log_depth_residualized")
meta <- fread(file.path(IN_DIR, "sample_metadata_stage1.tsv"))
depth_pca <- fread(file.path(IN_DIR, "input_variant_depth_pca_summary.tsv"))

summary_rows <- list()
overlap_rows <- list()
pres_rows <- list()
conc_rows <- list()

log_msg(sprintf("Input sensitivity settings: n_perm=%d", n_perm))

for (variant_name in variants) {
  mat <- readRDS(file.path(IN_DIR, paste0(variant_name, ".rds")))
  expr_by_core <- split_by_core(mat, meta)

  for (i in seq_len(nrow(settings))) {
    par <- settings[i]
    sid <- par$setting_id
    out_dir <- file.path(OUT, variant_name, sid)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    log_msg(sprintf("Evaluating %s / %s", variant_name, sid))

    net <- fit_consensus(expr_by_core, par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize)
    mods <- data.table(taxon = names(net$colors), module = as.character(net$colors))
    fwrite(mods, file.path(out_dir, "module_assignments.tsv"), sep = "\t")

    ref_path <- here("networkQC", "results", "full_eval", sid, "module_assignments.tsv")
    ref_mods <- fread(ref_path)
    comp <- compare_assignments(ref_mods, mods)
    ref_detail <- comp$ref_module_rows[[1]]
    ref_detail[, `:=`(variant = variant_name, setting_id = sid)]
    overlap_rows[[length(overlap_rows) + 1]] <- ref_detail

    pooled <- do.call(rbind, lapply(PARAMS$stage1_cores, function(core_id) expr_by_core[[core_id]]))
    MEs_train <- orderMEs(moduleEigengenes(pooled, net$colors)$eigengenes)
    MEs_valid <- orderMEs(moduleEigengenes(expr_by_core[[PARAMS$validation_core]], net$colors)$eigengenes)
    MEs_all <- rbind(
      as.data.table(MEs_train, keep.rownames = "sample"),
      as.data.table(MEs_valid, keep.rownames = "sample")
    )
    fwrite(MEs_all, file.path(out_dir, "module_eigengenes.tsv"), sep = "\t")

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
    fwrite(pres, file.path(out_dir, "preservation.tsv"), sep = "\t")
    pres_rows[[length(pres_rows) + 1]] <- pres

    conc <- age_aligned_concordance(MEs_all, meta)
    conc[, `:=`(variant = variant_name, setting_id = sid)]
    fwrite(conc, file.path(out_dir, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
    conc_rows[[length(conc_rows) + 1]] <- conc

    bio_pres <- pres[module_type == "biological"]
    bio_conc <- conc[!grepl("grey|gold", module)]
    dpc <- depth_pca[depth_pca$variant == variant_name]
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      variant = variant_name,
      setting_id = sid,
      power = par$power,
      deepSplit = par$deepSplit,
      mergeCutHeight = par$mergeCutHeight,
      minModuleSize = par$minModuleSize,
      non_grey_modules = mods[module != "grey", uniqueN(module)],
      grey_pct = mods[module == "grey", .N] / nrow(mods) * 100,
      mean_ref_best_jaccard = comp$mean_ref_best_jaccard,
      median_ref_best_jaccard = comp$median_ref_best_jaccard,
      min_ref_best_jaccard = comp$min_ref_best_jaccard,
      grey_delta_pct = comp$grey_delta_pct,
      same_grey_fraction = comp$same_grey_fraction,
      bio_pres_strong = bio_pres[preserved == "strong", .N],
      bio_pres_moderate = bio_pres[preserved == "moderate", .N],
      mean_concordance_pearson = mean(bio_conc$pearson_r, na.rm = TRUE),
      mean_concordance_spearman = mean(bio_conc$spearman_rho, na.rm = TRUE),
      mean_concordance_rmse = mean(bio_conc$rmse, na.rm = TRUE),
      cor_PC1_log_total = dpc$cor_PC1_log_total,
      cor_PC2_log_total = dpc$cor_PC2_log_total,
      cor_PC3_log_total = dpc$cor_PC3_log_total
    )
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
overlap_dt <- rbindlist(overlap_rows, fill = TRUE)
pres_dt <- rbindlist(pres_rows, fill = TRUE)
conc_dt <- rbindlist(conc_rows, fill = TRUE)

summary_dt[, drastic_change := (
  variant != "current_taxon_centered_log" &
    (mean_ref_best_jaccard < 0.50 |
       min_ref_best_jaccard < 0.20 |
       abs(grey_delta_pct) > 10 |
       bio_pres_strong < 4 |
       mean_concordance_pearson < 0.60)
)]
setorder(summary_dt, setting_id, variant)

fwrite(summary_dt, file.path(OUT, "input_sensitivity_summary.tsv"), sep = "\t")
fwrite(overlap_dt, file.path(OUT, "module_overlap_to_current.tsv"), sep = "\t")
fwrite(pres_dt, file.path(OUT, "preservation_all.tsv"), sep = "\t")
fwrite(conc_dt, file.path(OUT, "eigengene_concordance_all.tsv"), sep = "\t")

metric_cols <- c(
  "grey_pct", "mean_ref_best_jaccard", "min_ref_best_jaccard",
  "bio_pres_strong", "bio_pres_moderate",
  "mean_concordance_pearson", "mean_concordance_spearman", "mean_concordance_rmse",
  "cor_PC1_log_total"
)
heat <- melt(summary_dt, id.vars = c("variant", "setting_id"), measure.vars = metric_cols,
             variable.name = "metric", value.name = "value")
heat[, scaled := {
  r <- range(value, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) rep(0.5, .N) else (value - r[1]) / (r[2] - r[1])
}, by = metric]
heat[metric %in% c("grey_pct", "mean_concordance_rmse", "cor_PC1_log_total"), scaled := 1 - abs(scaled)]
heat[metric == "cor_PC1_log_total", scaled := {
  mx <- max(abs(value), na.rm = TRUE)
  if (!is.finite(mx) || mx == 0) rep(0.5, .N) else 1 - abs(value) / mx
}]
heat[, label := sprintf("%.2f", value)]
heat[, row_label := paste(setting_id, variant, sep = " | ")]
fwrite(heat, file.path(OUT, "input_sensitivity_heatmap_long.tsv"), sep = "\t")

png(file.path(FIG, "input_sensitivity_heatmap.png"), width = 2500, height = 1600, res = 220)
print(
  ggplot(heat, aes(metric, row_label, fill = scaled)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = 2) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey90", name = "Scaled") +
    labs(
      title = "Input preprocessing sensitivity for exp3, exp4, and opt5",
      subtitle = sprintf("Preservation permutations per run: %d; values printed in original units", n_perm),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

report <- file.path(BASE, "INPUT_EVALUATION_REBUILD_REPORT.md")
best_exp3 <- summary_dt[setting_id == "exp3"]
sink(report)
cat("# NetworkQC Input Evaluation Rebuild Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: focused preprocessing/depth sensitivity for `exp3`, `exp4`, and `opt5`\n")
cat("- Preservation permutations per sensitivity run: `", n_perm, "`\n\n", sep = "")

cat("## Input Variants\n\n")
cat("- `current_taxon_centered_log`: current main-pipeline WGCNA input, taxon-centered log raw counts.\n")
cat("- `sample_clr_raw`: sample-wise CLR on raw filtered counts with pseudocount 0.5.\n")
cat("- `deseq_length_log`: DESeq2 poscounts + reference-length normalized log counts, taxon-centered.\n")
cat("- `log_depth_residualized`: current log matrix after removing linear log-total-read signal per taxon.\n\n")

cat("## Summary Table\n\n")
cat("|setting|variant|grey_pct|mean_overlap|min_overlap|strong|moderate|pearson|spearman|rmse|PC1_depth|drastic|\n")
cat("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  cat(sprintf("|%s|%s|%.2f|%.3f|%.3f|%d|%d|%.3f|%.3f|%.3f|%.3f|%s|\n",
              r$setting_id, r$variant, r$grey_pct, r$mean_ref_best_jaccard, r$min_ref_best_jaccard,
              r$bio_pres_strong, r$bio_pres_moderate, r$mean_concordance_pearson,
              r$mean_concordance_spearman, r$mean_concordance_rmse, r$cor_PC1_log_total,
              ifelse(isTRUE(r$drastic_change), "yes", "no")))
}
cat("\n")

cat("## Initial Interpretation\n\n")
flagged <- summary_dt[drastic_change == TRUE]
if (nrow(flagged) == 0) {
  cat("No non-current variant crossed the pre-set drastic-change thresholds. This would support keeping `exp3` as the candidate network while moving to kME/topology checks.\n\n")
} else {
  cat("At least one non-current variant crossed the pre-set drastic-change thresholds. Review these before locking the current network:\n\n")
  for (i in seq_len(nrow(flagged))) {
    r <- flagged[i]
    cat(sprintf("- `%s / %s`: mean overlap %.3f, min overlap %.3f, grey delta %.2f, strong biological modules %d, Pearson %.3f\n",
                r$setting_id, r$variant, r$mean_ref_best_jaccard, r$min_ref_best_jaccard,
                r$grey_delta_pct, r$bio_pres_strong, r$mean_concordance_pearson))
  }
  cat("\n")
}

cat("For `exp3`, compare the non-current rows against `current_taxon_centered_log`. If module overlap stays high and preservation/concordance remain strong, the current `exp3` decision is robust to input preprocessing. If overlap collapses or depth-controlled variants change the module structure, rerun the broader parameter sweep on the preferred corrected input.\n\n")

cat("## Output Map\n\n")
cat("- Input variants: `networkQC/input_evaluation/results/inputs/*.rds`\n")
cat("- Summary table: `networkQC/input_evaluation/results/sensitivity/input_sensitivity_summary.tsv`\n")
cat("- Module overlap details: `networkQC/input_evaluation/results/sensitivity/module_overlap_to_current.tsv`\n")
cat("- Preservation details: `networkQC/input_evaluation/results/sensitivity/preservation_all.tsv`\n")
cat("- Concordance details: `networkQC/input_evaluation/results/sensitivity/eigengene_concordance_all.tsv`\n")
cat("- Heatmap: `networkQC/input_evaluation/results/figures/input_sensitivity_heatmap.png`\n")
sink()

log_msg("Input sensitivity complete. Report: ", report)
