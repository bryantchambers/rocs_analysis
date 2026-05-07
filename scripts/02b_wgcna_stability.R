#!/usr/bin/env Rscript
# 02b_wgcna_stability.R — WGCNA stability diagnostics via bootstrap consensus reruns
#
# Purpose:
#   1) quantify module assignment stability under resampling
#   2) summarize preservation + eigengene concordance outputs from 02_wgcna.R
#   3) write a novice-readable markdown report
#
# Outputs: results/wgcna_stability/
#   module_stability_bootstrap.tsv
#   module_stability_summary.tsv
#   module_size_sensitivity.tsv
#   WGCNA_STABILITY_REPORT.md

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

args <- commandArgs(trailingOnly = TRUE)
run_mode <- PARAMS$wgcna_run_mode
if (length(args) > 0) {
  mode_arg <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
  if (!is.na(mode_arg) && nzchar(mode_arg)) run_mode <- mode_arg
}
if (!run_mode %in% c("build", "final")) {
  stop("Invalid WGCNA mode: ", run_mode, ". Use --mode=build or --mode=final.")
}

n_boot <- if (run_mode == "final") PARAMS$wgcna_stability_bootstrap_final else PARAMS$wgcna_stability_bootstrap_build
WOUT <- file.path(RESULTS$stage1, "wgcna")
SOUT <- RESULTS$wgcna_stability
dir.create(SOUT, recursive = TRUE, showWarnings = FALSE)

log_msg(sprintf("Stability mode: %s, bootstrap replicates: %d", run_mode, n_boot))

# ---- load baseline objects/results -------------------------------------------
vst  <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
mods <- fread(file.path(WOUT, "module_assignments.tsv"))
soft <- fread(file.path(WOUT, "soft_power.tsv"))
pres <- fread(file.path(WOUT, "preservation.tsv"))
pres_bio <- fread(file.path(WOUT, "preservation_biological.tsv"))
conc <- fread(file.path(WOUT, "eigengene_concordance_age_aligned.tsv"))
setnames(conc, old = names(conc), new = sub("\\.V1$", "", names(conc)))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- meta[core == core_id, label]
  samps <- intersect(samps, rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

# Recover selected soft power from step-02 rule.
soft[, signedR2 := ifelse(slope < 0, SFT.R.sq, 0)]
passing <- soft[signedR2 >= 0.80, .(min_pass = min(Power)), by = core]
if (nrow(passing) == length(PARAMS$stage1_cores)) {
  soft_power <- max(passing$min_pass)
} else {
  best <- soft[slope < 0 & mean.k. > 5, .(best_r2 = max(signedR2)), by = Power]
  soft_power <- if (nrow(best) > 0) best[which.max(best_r2), Power] else 12L
}

baseline_assign <- setNames(mods$module, mods$taxon)
bio_modules <- setdiff(sort(unique(mods$module)), c("grey", "gold"))

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(NA_real_)
  inter / union
}

bootstrap_once <- function(seed_offset) {
  set.seed(PARAMS$seed + seed_offset)
  sampled <- lapply(PARAMS$stage1_cores, function(core_id) {
    dat <- expr_by_core[[core_id]]
    idx <- sample(seq_len(nrow(dat)), size = nrow(dat), replace = TRUE)
    dat[idx, , drop = FALSE]
  })
  names(sampled) <- PARAMS$stage1_cores

  multiExpr <- lapply(PARAMS$stage1_cores, function(core_id) list(data = sampled[[core_id]]))
  names(multiExpr) <- PARAMS$stage1_cores

  fit <- tryCatch(
    blockwiseConsensusModules(
      multiExpr      = multiExpr,
      power          = soft_power,
      networkType    = "signed",
      corType        = "pearson",
      maxBlockSize   = 5000,
      minModuleSize  = PARAMS$wgcna_min_module_size,
      deepSplit      = PARAMS$wgcna_deep_split,
      mergeCutHeight = PARAMS$wgcna_merge_cut_height,
      numericLabels  = FALSE,
      saveTOMs       = FALSE,
      verbose        = 0
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  asg <- fit$colors
  names(asg) <- colnames(sampled[[1]])
  asg
}

log_msg("Running bootstrap stability reruns...")
boot_assign <- vector("list", n_boot)
for (i in seq_len(n_boot)) {
  if (i %% 10 == 0 || i == 1 || i == n_boot) log_msg(sprintf("  replicate %d/%d", i, n_boot))
  boot_assign[[i]] <- bootstrap_once(i)
}
ok <- which(vapply(boot_assign, Negate(is.null), logical(1)))
if (length(ok) == 0) stop("All bootstrap reruns failed; no stability output produced.")
boot_assign <- boot_assign[ok]
log_msg(sprintf("Bootstrap successful runs: %d/%d", length(boot_assign), n_boot))

stability_rows <- rbindlist(lapply(seq_along(boot_assign), function(i) {
  asg <- boot_assign[[i]]
  mods_here <- setNames(asg[names(baseline_assign)], names(baseline_assign))
  boot_mod_names <- setdiff(sort(unique(mods_here)), c("grey", "gold"))
  rbindlist(lapply(bio_modules, function(m) {
    base_genes <- names(baseline_assign)[baseline_assign == m]
    if (length(boot_mod_names) == 0) {
      return(data.table(
        replicate = i,
        module = m,
        matched_boot_module = NA_character_,
        jaccard = NA_real_,
        baseline_size = length(base_genes),
        bootstrap_size = NA_integer_
      ))
    }
    j_by_boot <- rbindlist(lapply(boot_mod_names, function(bm) {
      bg <- names(mods_here)[mods_here == bm]
      data.table(boot_module = bm, jaccard = jaccard(base_genes, bg), bootstrap_size = length(bg))
    }))
    best <- j_by_boot[which.max(jaccard)]
    data.table(
      replicate = i,
      module = m,
      matched_boot_module = best$boot_module,
      jaccard = best$jaccard,
      baseline_size = length(base_genes),
      bootstrap_size = best$bootstrap_size
    )
  }))
}), fill = TRUE)

stability_sum <- stability_rows[, .(
  n = .N,
  jaccard_median = median(jaccard, na.rm = TRUE),
  jaccard_p05 = quantile(jaccard, 0.05, na.rm = TRUE),
  jaccard_p95 = quantile(jaccard, 0.95, na.rm = TRUE),
  bootstrap_size_median = median(bootstrap_size, na.rm = TRUE)
), by = module][order(-jaccard_median)]

size_sens <- stability_rows[, .(
  size_ratio_median = median(bootstrap_size / baseline_size, na.rm = TRUE),
  size_ratio_p05 = quantile(bootstrap_size / baseline_size, 0.05, na.rm = TRUE),
  size_ratio_p95 = quantile(bootstrap_size / baseline_size, 0.95, na.rm = TRUE)
), by = module][order(module)]

fwrite(stability_rows, file.path(SOUT, "module_stability_bootstrap.tsv"), sep = "\t")
fwrite(stability_sum, file.path(SOUT, "module_stability_summary.tsv"), sep = "\t")
fwrite(size_sens, file.path(SOUT, "module_size_sensitivity.tsv"), sep = "\t")

# ---- novice-readable markdown report -----------------------------------------
report <- file.path(SOUT, "WGCNA_STABILITY_REPORT.md")
sink(report)
cat("# WGCNA Stability Report\n\n")
cat("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Mode: `", run_mode, "`\n", sep = "")
cat("- Bootstrap reruns requested: ", n_boot, "\n", sep = "")
cat("- Bootstrap reruns successful: ", length(boot_assign), "\n\n", sep = "")

cat("## What This Means\n\n")
cat("This report checks whether WGCNA module assignments stay similar when samples are re-drawn.\n")
cat("Higher Jaccard overlap means module membership is more stable.\n")
cat("We also summarize holdout preservation and age-aligned eigengene concordance from `02_wgcna.R`.\n\n")

cat("## 1) Module Stability Under Bootstrap Resampling\n\n")
cat("|module|median_jaccard|p05|p95|baseline_size|bootstrap_size_median|\n")
cat("|---|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(stability_sum))) {
  m <- stability_sum$module[i]
  base_size <- unique(stability_rows[module == m, baseline_size])[1]
  cat(sprintf("|%s|%.3f|%.3f|%.3f|%d|%.1f|\n",
              m, stability_sum$jaccard_median[i], stability_sum$jaccard_p05[i],
              stability_sum$jaccard_p95[i], base_size, stability_sum$bootstrap_size_median[i]))
}
cat("\n")

cat("Interpretation guide:\n")
cat("- `median_jaccard >= 0.70`: strong stability\n")
cat("- `0.50 to 0.70`: moderate stability\n")
cat("- `< 0.50`: potentially unstable module definition\n\n")

cat("## 2) Preservation Results (R1 -> R2)\n\n")
cat("All modules (`preservation.tsv`):\n\n")
cat("|module|Zsummary|preserved|module_type|\n|---|---:|---|---|\n")
for (i in seq_len(nrow(pres))) {
  mt <- if ("module_type" %in% names(pres)) pres$module_type[i] else ifelse(pres$module[i] %in% c("grey", "gold"), "technical", "biological")
  cat(sprintf("|%s|%.3f|%s|%s|\n", pres$module[i], pres$Zsummary[i], pres$preserved[i], mt))
}
cat("\nBiological-only modules (`preservation_biological.tsv`):\n\n")
cat("|module|Zsummary|preserved|\n|---|---:|---|\n")
for (i in seq_len(nrow(pres_bio))) {
  cat(sprintf("|%s|%.3f|%s|\n", pres_bio$module[i], pres_bio$Zsummary[i], pres_bio$preserved[i]))
}
cat("\n")

cat("## 3) Age-Aligned Eigengene Concordance (R1 vs R2)\n\n")
cat("Computed on a common age grid after interpolation.\n\n")
cat("|module|pearson_r|spearman_rho|rmse|\n|---|---:|---:|---:|\n")
for (i in seq_len(nrow(conc))) {
  cat(sprintf("|%s|%.3f|%.3f|%.3f|\n", conc$module[i], conc$pearson_r[i], conc$spearman_rho[i], conc$rmse[i]))
}
cat("\n")

cat("Interpretation guide:\n")
cat("- Higher positive `pearson_r` / `spearman_rho` indicates better cross-core agreement.\n")
cat("- Lower `rmse` indicates closer eigengene trajectories.\n\n")

cat("## 4) Recommended Next Action\n\n")
cat("Use `--mode=build` during daytime development, and `--mode=final` overnight.\n")
cat("Treat modules with low bootstrap Jaccard or low R1/R2 concordance as lower-confidence for downstream biological claims.\n")
sink()

log_msg("Saved stability outputs to ", SOUT)
log_msg("Saved report: ", report)
