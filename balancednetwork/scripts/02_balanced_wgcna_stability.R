#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
})

source(here("config.R"))
allowWGCNAThreads()
set.seed(PARAMS$seed)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA_character_) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}

run_mode <- arg_val("mode", "build")
if (!run_mode %in% c("build", "final")) stop("Invalid mode.")
n_boot <- as.integer(arg_val("n_boot", ifelse(run_mode == "final", "1000", "100")))
force <- as.integer(arg_val("force", "0")) == 1L
if (!is.finite(n_boot) || n_boot <= 0) stop("Invalid n_boot.")

OUT_BASE <- here("balancednetwork", "results")
OUT_TAB <- file.path(OUT_BASE, "tables")
OUT_WGCNA <- file.path(OUT_BASE, "wgcna")
OUT_STAB <- file.path(OUT_BASE, "stability")
OUT_LOG <- file.path(OUT_BASE, "logs")
dir.create(OUT_STAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOG, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(OUT_STAB, "module_stability_summary.tsv")
if (file.exists(summary_file) && !force) {
  message("[balancednetwork] stability exists; skipping (use --force=1 to rebuild).")
  quit(save = "no")
}

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
mods <- fread(file.path(OUT_WGCNA, "module_assignments.tsv"))
sel <- fread(file.path(OUT_TAB, "balanced_baseline_samples.tsv"))
bin_quota <- fread(file.path(OUT_TAB, "balance_bin_quotas.tsv"))

train_cores <- PARAMS$stage1_cores
design_pool <- merge(
  sel[, .(label, core, age_bin)],
  sel[, .(label)],
  by = "label"
)
design_pool <- unique(design_pool)

pool_by_core_bin <- lapply(train_cores, function(core_id) {
  dt <- sel[core == core_id, .(label, age_bin)]
  split(dt$label, dt$age_bin)
})
names(pool_by_core_bin) <- train_cores

baseline_assign <- setNames(mods$module, mods$taxon)
bio_modules <- setdiff(sort(unique(mods$module)), c("grey", "gold"))

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

sample_manifest_rows <- vector("list", n_boot)
boot_rows <- vector("list", n_boot)
fail_rows <- list()

for (i in seq_len(n_boot)) {
  if (i %% 25 == 0 || i == 1 || i == n_boot) {
    message(sprintf("[balancednetwork] bootstrap %d/%d", i, n_boot))
  }
  set.seed(PARAMS$seed + i)

  sampled_labels <- lapply(train_cores, function(core_id) {
    core_pools <- pool_by_core_bin[[core_id]]
    taken <- unlist(lapply(seq_len(nrow(bin_quota)), function(j) {
      bin_id <- bin_quota$age_bin[j]
      quota <- bin_quota$quota_per_core[j]
      pool <- core_pools[[bin_id]]
      if (is.null(pool) || length(pool) == 0 || quota <= 0) return(character(0))
      sample(pool, size = quota, replace = TRUE)
    }), use.names = FALSE)
    taken
  })
  names(sampled_labels) <- train_cores

  sample_manifest_rows[[i]] <- rbindlist(lapply(train_cores, function(core_id) {
    data.table(replicate = i, core = core_id, label = sampled_labels[[core_id]])
  }))

  multiExpr <- lapply(train_cores, function(core_id) {
    list(data = vst[sampled_labels[[core_id]], names(baseline_assign), drop = FALSE])
  })
  names(multiExpr) <- train_cores

  fit <- tryCatch(
    blockwiseConsensusModules(
      multiExpr = multiExpr,
      power = 12,
      networkType = "signed",
      corType = "pearson",
      maxBlockSize = 5000,
      minModuleSize = 20,
      deepSplit = 3,
      mergeCutHeight = 0.25,
      numericLabels = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    fail_rows[[length(fail_rows) + 1L]] <- data.table(replicate = i, error = conditionMessage(fit))
    next
  }

  asg <- setNames(as.character(fit$colors), names(fit$colors))
  asg <- asg[names(baseline_assign)]
  boot_mods <- setdiff(sort(unique(asg)), c("grey", "gold"))

  boot_rows[[i]] <- rbindlist(lapply(bio_modules, function(m) {
    bg <- names(baseline_assign)[baseline_assign == m]
    if (!length(boot_mods)) {
      return(data.table(replicate = i, module = m, matched_boot_module = NA_character_, jaccard = NA_real_))
    }
    js <- sapply(boot_mods, function(bm) jaccard(bg, names(asg)[asg == bm]))
    best_idx <- which.max(js)
    data.table(replicate = i, module = m, matched_boot_module = boot_mods[best_idx], jaccard = js[best_idx])
  }))
}

sample_manifest <- rbindlist(sample_manifest_rows, fill = TRUE)
boot_dt <- rbindlist(boot_rows, fill = TRUE)
fail_dt <- if (length(fail_rows)) rbindlist(fail_rows) else data.table(replicate = integer(), error = character())

if (nrow(boot_dt) == 0) stop("All balanced bootstrap replicates failed.")

stab_summary <- boot_dt[, .(
  n = .N,
  jaccard_median = median(jaccard, na.rm = TRUE),
  jaccard_p05 = quantile(jaccard, 0.05, na.rm = TRUE),
  jaccard_p95 = quantile(jaccard, 0.95, na.rm = TRUE)
), by = module][order(-jaccard_median)]

run_summary <- data.table(
  metric = c("mode", "n_boot_requested", "n_boot_successful", "n_boot_failed"),
  value = c(run_mode, as.character(n_boot), as.character(length(unique(boot_dt$replicate))), as.character(nrow(fail_dt)))
)

orig_assign_file <- file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv")
if (file.exists(orig_assign_file)) {
  orig <- fread(orig_assign_file)
  orig_asg <- setNames(orig$module, orig$taxon)
  common <- intersect(names(orig_asg), names(baseline_assign))
  common <- common[is.finite(match(common, names(baseline_assign)))]
  orig_mods <- setdiff(sort(unique(orig_asg[common])), c("grey", "gold"))
  bal_mods <- setdiff(sort(unique(baseline_assign[common])), c("grey", "gold"))
  overlap <- rbindlist(lapply(bal_mods, function(m) {
    bg <- names(baseline_assign)[baseline_assign == m]
    bg <- intersect(bg, common)
    if (!length(orig_mods)) return(data.table(balanced_module = m, best_original_module = NA_character_, best_jaccard = NA_real_))
    js <- sapply(orig_mods, function(om) jaccard(bg, names(orig_asg)[orig_asg == om]))
    idx <- which.max(js)
    data.table(balanced_module = m, best_original_module = orig_mods[idx], best_jaccard = js[idx])
  }))
  fwrite(overlap, file.path(OUT_STAB, "balanced_vs_original_module_overlap.tsv"), sep = "\t")
}

fwrite(sample_manifest, file.path(OUT_STAB, "bootstrap_sample_manifest.tsv"), sep = "\t")
fwrite(boot_dt, file.path(OUT_STAB, "module_stability_bootstrap.tsv"), sep = "\t")
fwrite(stab_summary, summary_file, sep = "\t")
fwrite(fail_dt, file.path(OUT_STAB, "bootstrap_failures.tsv"), sep = "\t")
fwrite(run_summary, file.path(OUT_STAB, "stability_run_summary.tsv"), sep = "\t")

message(sprintf(
  "[balancednetwork] stability complete: requested=%d successful=%d failed=%d",
  n_boot, length(unique(boot_dt$replicate)), nrow(fail_dt)
))
