#!/usr/bin/env Rscript
# 02b_wgcna_stability.R - WGCNA module stability diagnostics for the M/Bryant merge
#
# This stage leaves module construction untouched and evaluates whether the
# selected modules are stable under sample resampling. It is intentionally placed
# after 02_wgcna_main.R and before M's HMM stage.

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_mode <- PARAMS$wgcna_run_mode
if (length(args) > 0) {
  mode_arg <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
  if (!is.na(mode_arg) && nzchar(mode_arg)) run_mode <- mode_arg
}
if (!run_mode %in% c("build", "final")) {
  stop("Invalid WGCNA mode: ", run_mode, ". Use --mode=build or --mode=final.")
}

n_boot <- if (run_mode == "final") {
  PARAMS$wgcna_stability_bootstrap_final
} else {
  PARAMS$wgcna_stability_bootstrap_build
}
ensure_dirs(c(DIRS$wgcna_stability))
log_msg(sprintf("WGCNA stability mode: %s; bootstrap replicates: %d", run_mode, n_boot))

clr <- readRDS(file.path(DIRS$main, "clr_matrix_train_centered.rds"))
meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
mods <- fread(file.path(DIRS$main, "module_assignments.tsv"))
soft <- fread(file.path(DIRS$main, "soft_power_scan.tsv"))
pres <- fread(file.path(DIRS$main, "module_preservation_validation.tsv"))
pres_bio_path <- file.path(DIRS$main, "module_preservation_validation_biological.tsv")
pres_bio <- if (file.exists(pres_bio_path)) fread(pres_bio_path) else pres[!module %in% c("grey", "gold")]
conc <- fread(file.path(DIRS$main, "eigengene_concordance_age_aligned.tsv"))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  ids <- meta_main[core == core_id, label]
  ids <- intersect(ids, rownames(clr))
  clr[ids, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

training_manifest_path <- file.path(DIRS$main, "wgcna_training_samples.tsv")
if (!file.exists(training_manifest_path)) {
  stop("Missing WGCNA training manifest: ", training_manifest_path, ". Run 01_data_prep.R first.")
}
wgcna_training_manifest <- fread(training_manifest_path)
fit_expr_by_core <- lapply(PARAMS$training_cores, function(core_id) {
  ids <- wgcna_training_manifest[core == core_id, label]
  ids <- intersect(ids, rownames(clr))
  if (length(ids) == 0) stop("No selected WGCNA training samples for core: ", core_id)
  clr[ids, , drop = FALSE]
})
names(fit_expr_by_core) <- PARAMS$training_cores

balance_quotas <- data.table()
pool_by_core_bin <- list()
if (identical(PARAMS$wgcna_input_strategy, "balanced")) {
  quota_path <- file.path(DIRS$main, "balance_bin_quotas.tsv")
  if (!file.exists(quota_path)) stop("Missing balanced quota table: ", quota_path)
  balance_quotas <- fread(quota_path)
  if (!all(c("age_bin", "quota_per_core") %in% names(balance_quotas))) {
    stop("Balanced quota table must contain age_bin and quota_per_core.")
  }
  if (!"age_bin" %in% names(wgcna_training_manifest)) {
    stop("Balanced WGCNA training manifest must contain age_bin.")
  }
  pool_by_core_bin <- lapply(PARAMS$training_cores, function(core_id) {
    dt <- wgcna_training_manifest[core == core_id, .(label, age_bin)]
    split(dt$label, dt$age_bin)
  })
  names(pool_by_core_bin) <- PARAMS$training_cores
}

if (!is.null(PARAMS$wgcna_soft_power) && is.finite(PARAMS$wgcna_soft_power)) {
  soft_power <- as.integer(PARAMS$wgcna_soft_power)
} else {
  soft[, signedR2 := ifelse(slope < 0, SFT.R.sq, 0)]
  passing <- soft[signedR2 >= PARAMS$soft_power_target_r2, .(min_pass = min(Power)), by = core]
  if (nrow(passing) == length(PARAMS$training_cores)) {
    soft_power <- max(passing$min_pass)
  } else {
    best <- soft[slope < 0 & mean.k. > 5, .(best_r2 = max(signedR2)), by = Power]
    soft_power <- if (nrow(best) > 0) best[which.max(best_r2), Power] else 12L
  }
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
  if (identical(PARAMS$wgcna_input_strategy, "balanced")) {
    sampled <- lapply(PARAMS$training_cores, function(core_id) {
      core_pools <- pool_by_core_bin[[core_id]]
      sampled_labels <- unlist(lapply(seq_len(nrow(balance_quotas)), function(j) {
        bin_id <- balance_quotas$age_bin[j]
        quota <- balance_quotas$quota_per_core[j]
        pool <- core_pools[[bin_id]]
        if (is.null(pool) || length(pool) == 0 || quota <= 0) return(character(0))
        sample(pool, size = quota, replace = TRUE)
      }), use.names = FALSE)
      clr[sampled_labels, names(baseline_assign), drop = FALSE]
    })
  } else {
    sampled <- lapply(PARAMS$training_cores, function(core_id) {
      dat <- fit_expr_by_core[[core_id]]
      idx <- sample(seq_len(nrow(dat)), size = nrow(dat), replace = TRUE)
      dat[idx, , drop = FALSE]
    })
  }
  names(sampled) <- PARAMS$training_cores

  multiExpr <- lapply(PARAMS$training_cores, function(core_id) list(data = sampled[[core_id]]))
  names(multiExpr) <- PARAMS$training_cores

  fit <- tryCatch(
    blockwiseConsensusModules(
      multiExpr = multiExpr,
      power = soft_power,
      networkType = "signed",
      corType = "pearson",
      maxBlockSize = 5000,
      minModuleSize = PARAMS$wgcna_min_module_size,
      deepSplit = PARAMS$wgcna_deep_split,
      mergeCutHeight = PARAMS$wgcna_merge_cut_height,
      numericLabels = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  asg <- fit$colors
  names(asg) <- colnames(sampled[[1]])
  asg
}

log_msg("Running bootstrap WGCNA stability reruns...")
boot_assign <- vector("list", n_boot)
for (i in seq_len(n_boot)) {
  if (i == 1 || i == n_boot || i %% 10 == 0) {
    log_msg(sprintf("Bootstrap replicate %d/%d", i, n_boot))
  }
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

fwrite(stability_rows, file.path(DIRS$wgcna_stability, "module_stability_bootstrap.tsv"), sep = "\t")
fwrite(stability_sum, file.path(DIRS$wgcna_stability, "module_stability_summary.tsv"), sep = "\t")
fwrite(size_sens, file.path(DIRS$wgcna_stability, "module_size_sensitivity.tsv"), sep = "\t")

report <- file.path(DIRS$wgcna_stability, "WGCNA_STABILITY_REPORT.md")
sink(report)
cat("# WGCNA Stability Report\n\n")
cat("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Mode: `", run_mode, "`\n", sep = "")
cat("- Bootstrap reruns requested: ", n_boot, "\n", sep = "")
cat("- Bootstrap reruns successful: ", length(boot_assign), "\n", sep = "")
cat("- Soft power: ", soft_power, "\n", sep = "")
cat("- deepSplit: ", PARAMS$wgcna_deep_split, "\n", sep = "")
cat("- mergeCutHeight: ", PARAMS$wgcna_merge_cut_height, "\n", sep = "")
cat("- minModuleSize: ", PARAMS$wgcna_min_module_size, "\n\n", sep = "")
cat("- WGCNA input strategy: `", PARAMS$wgcna_input_strategy, "`\n", sep = "")
cat("- WGCNA profile: `", PARAMS$wgcna_profile, "`\n", sep = "")
cat("- HMM input balancing status: `", PARAMS$hmm_input_balancing_status, "`\n\n", sep = "")

cat("## What This Checks\n\n")
cat("This report asks whether taxa stay in similar biological modules when the training samples are resampled with replacement.\n")
cat("Higher Jaccard overlap means more reproducible module membership. Preservation and age-aligned R1/R2 concordance are summarized from `02_wgcna_main.R`.\n\n")

cat("## 1) Bootstrap Module Stability\n\n")
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

cat("## 2) Preservation Results (GeoB25202_R1 -> GeoB25202_R2)\n\n")
cat("All modules:\n\n")
cat("|module|Zsummary|Zdensity|Zconnectivity|preserved|module_type|\n")
cat("|---|---:|---:|---:|---|---|\n")
for (i in seq_len(nrow(pres))) {
  mt <- if ("module_type" %in% names(pres)) pres$module_type[i] else ifelse(pres$module[i] %in% c("grey", "gold"), "technical", "biological")
  cat(sprintf("|%s|%.3f|%.3f|%.3f|%s|%s|\n",
              pres$module[i], pres$Zsummary[i], pres$Zdensity[i], pres$Zconnectivity[i], pres$preserved[i], mt))
}
cat("\nBiological modules only:\n\n")
cat("|module|Zsummary|preserved|\n|---|---:|---|\n")
for (i in seq_len(nrow(pres_bio))) {
  cat(sprintf("|%s|%.3f|%s|\n", pres_bio$module[i], pres_bio$Zsummary[i], pres_bio$preserved[i]))
}
cat("\n")

cat("## 3) Age-Aligned Eigengene Concordance\n\n")
cat("R1 and R2 are interpolated onto a common age grid before correlation.\n\n")
cat("|module|pearson_r|spearman_rho|rmse|\n|---|---:|---:|---:|\n")
for (i in seq_len(nrow(conc))) {
  cat(sprintf("|%s|%.3f|%.3f|%.3f|\n", conc$module[i], conc$pearson_r[i], conc$spearman_rho[i], conc$rmse[i]))
}
cat("\n")

cat("## Practical Interpretation\n\n")
cat("- Use `--mode=build` for development checks.\n")
cat("- Use `--mode=final` for final claims or overnight runs.\n")
cat("- Treat low-Jaccard modules as lower confidence in downstream biological interpretation.\n")
sink()

write_run_metadata(
  file.path(DIRS$wgcna_stability, "02b_wgcna_stability_run_metadata.tsv"),
  "02b_wgcna_stability.R",
  extra = list(
    run_mode = run_mode,
    bootstrap_replicates_requested = n_boot,
    bootstrap_replicates_successful = length(boot_assign),
    wgcna_input_strategy = PARAMS$wgcna_input_strategy,
    wgcna_profile = PARAMS$wgcna_profile,
    hmm_input_balancing_status = PARAMS$hmm_input_balancing_status,
    soft_power = soft_power,
    min_module_size = PARAMS$wgcna_min_module_size,
    deep_split = PARAMS$wgcna_deep_split,
    merge_cut_height = PARAMS$wgcna_merge_cut_height
  )
)
write_session_info(file.path(DIRS$wgcna_stability, "02b_wgcna_stability_sessionInfo.txt"))

log_msg("Saved WGCNA stability outputs to ", DIRS$wgcna_stability)
