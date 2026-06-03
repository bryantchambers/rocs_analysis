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
arg_val <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}

n_boot <- as.integer(arg_val("n_boot", as.character(PARAMS$wgcna_stability_bootstrap_build)))
n_perm <- as.integer(arg_val("n_perm", as.character(PARAMS$wgcna_preservation_permutations_build)))

BASE <- here("InputQC", "st8_low_taxa_network_sensitivity")
OUT_TABLE <- file.path(BASE, "results", "tables")
OUT_FIG <- file.path(BASE, "results", "figures")
OUT_SETTINGS <- file.path(BASE, "results", "settings")
REPORT <- file.path(BASE, "ST8_LOW_TAXA_NETWORK_SENSITIVITY_REPORT.md")
LOG <- file.path(BASE, "results", "st8_low_taxa_network_sensitivity.log")
PROGRESS <- file.path(BASE, "results", "progress.tsv")

dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SETTINGS, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) invisible(file.remove(LOG))
if (file.exists(PROGRESS)) invisible(file.remove(PROGRESS))

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}

progress_update <- function(setting_id, phase, status, details = "") {
  row <- data.table(
    ts = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    setting_id = setting_id,
    phase = phase,
    status = status,
    details = details
  )
  fwrite(row, PROGRESS, sep = "\t", append = file.exists(PROGRESS), col.names = !file.exists(PROGRESS))
}

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

rank_settings <- function(res) {
  score <- copy(res[is.finite(grey_pct) & is.finite(mean_bootstrap_jaccard)])
  norm <- function(x, high_better = TRUE) {
    r <- range(x, na.rm = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0.5, length(x)))
    z <- (x - r[1]) / (r[2] - r[1])
    if (!high_better) z <- 1 - z
    z
  }
  score[, s_grey := norm(grey_pct, FALSE)]
  score[, s_boot := norm(mean_bootstrap_jaccard, TRUE)]
  score[, s_pres := norm(bio_pres_strong + 0.5 * bio_pres_moderate, TRUE)]
  score[, s_conc := norm(mean_concordance_pearson, TRUE)]
  score[, s_balance := norm(mean_balanced_jaccard, TRUE)]
  score[, final_score := 0.30 * s_boot + 0.25 * s_pres + 0.20 * s_grey + 0.15 * s_conc + 0.10 * s_balance]
  setorder(score, -final_score)
  score[, sensitivity_rank := .I]
  score
}

age_aligned_concordance <- function(MEs_all, meta_dt, me_cols) {
  meta_age <- meta_dt[, .(label, core, age_kyr = y_bp / 1000)]
  dt <- merge(MEs_all, meta_age, by.x = "sample", by.y = "label", all.x = TRUE)
  r1 <- PARAMS$stage1_cores[[length(PARAMS$stage1_cores)]]
  r2 <- PARAMS$validation_core
  rbindlist(lapply(me_cols, function(me) {
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
      pearson_r = unname(cor(y1, y2, method = "pearson")),
      spearman_rho = unname(cor(y1, y2, method = "spearman")),
      rmse = sqrt(mean((y1 - y2)^2))
    )
  }))
}

log_msg(sprintf("Starting ST8 low-taxa removal network sensitivity: n_boot=%d, n_perm=%d", n_boot, n_perm))
vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
remove_samples <- fread(here("InputQC", "st8_low_taxa_review", "results", "tables", "st8_low_taxa_sample_table.tsv"))$sample
remove_samples <- intersect(remove_samples, rownames(vst))
keep_samples <- setdiff(rownames(vst), remove_samples)
vst_filt <- vst[keep_samples, , drop = FALSE]
meta_filt <- meta[label %in% keep_samples]

sample_summary <- rbindlist(list(
  meta[, .(n = .N), by = core][, dataset := "original"],
  meta_filt[, .(n = .N), by = core][, dataset := "filtered"]
), fill = TRUE)
fwrite(sample_summary, file.path(OUT_TABLE, "sample_counts_before_after.tsv"), sep = "\t")
fwrite(data.table(sample = remove_samples), file.path(OUT_TABLE, "removed_st8_low_taxa_samples.tsv"), sep = "\t")
log_msg(sprintf("Removed %d samples: %s", length(remove_samples), paste(remove_samples, collapse = ", ")))

settings <- data.table(
  setting_id = c("exp3", "exp4", "opt5"),
  power = c(12L, 12L, 12L),
  deepSplit = c(3L, 3L, 1L),
  mergeCutHeight = c(0.25, 0.20, 0.20),
  minModuleSize = c(20L, 20L, 20L)
)
fwrite(settings, file.path(OUT_TABLE, "candidate_settings.tsv"), sep = "\t")

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- intersect(meta_filt[core == core_id, label], rownames(vst_filt))
  vst_filt[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

fit_consensus <- function(power, deepSplit, mergeCutHeight, minModuleSize, sample_override = NULL) {
  mx <- lapply(PARAMS$stage1_cores, function(core_id) {
    dat <- expr_by_core[[core_id]]
    if (!is.null(sample_override)) dat <- dat[sample_override[[core_id]], , drop = FALSE]
    list(data = dat)
  })
  names(mx) <- PARAMS$stage1_cores
  blockwiseConsensusModules(
    multiExpr = mx,
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

evaluate_setting <- function(par) {
  sid <- par$setting_id
  sdir <- file.path(OUT_SETTINGS, sid)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  log_msg(sprintf("Evaluating %s", sid))
  progress_update(sid, "fit_consensus", "start")
  net <- fit_consensus(par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize)
  progress_update(sid, "fit_consensus", "ok")

  mods <- data.table(taxon = names(net$colors), module = as.character(net$colors))
  fwrite(mods, file.path(sdir, "module_assignments.tsv"), sep = "\t")

  progress_update(sid, "eigengenes", "start")
  pooled <- do.call(rbind, lapply(PARAMS$stage1_cores, function(c) expr_by_core[[c]]))
  MEs_train <- orderMEs(moduleEigengenes(pooled, net$colors)$eigengenes)
  MEs_valid <- orderMEs(moduleEigengenes(expr_by_core[[PARAMS$validation_core]], net$colors)$eigengenes)
  MEs_all <- rbind(
    as.data.table(MEs_train, keep.rownames = "sample"),
    as.data.table(MEs_valid, keep.rownames = "sample")
  )
  fwrite(MEs_all, file.path(sdir, "module_eigengenes.tsv"), sep = "\t")
  progress_update(sid, "eigengenes", "ok")

  progress_update(sid, "preservation", "start", sprintf("n_perm=%d", n_perm))
  r1 <- PARAMS$stage1_cores[[length(PARAMS$stage1_cores)]]
  r2 <- PARAMS$validation_core
  var_r1 <- apply(expr_by_core[[r1]], 2, var)
  var_r2 <- apply(expr_by_core[[r2]], 2, var)
  good <- names(which(var_r1 > 0 & var_r2 > 0))
  mp <- modulePreservation(
    multiData = list(
      R1 = list(data = expr_by_core[[r1]][, good]),
      R2 = list(data = expr_by_core[[r2]][, good])
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
  fwrite(pres, file.path(sdir, "preservation.tsv"), sep = "\t")
  fwrite(pres[module_type == "biological"], file.path(sdir, "preservation_biological.tsv"), sep = "\t")
  progress_update(sid, "preservation", "ok")

  progress_update(sid, "concordance", "start")
  conc <- age_aligned_concordance(MEs_all, meta_filt, setdiff(names(MEs_all), "sample"))
  fwrite(conc, file.path(sdir, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
  progress_update(sid, "concordance", "ok")

  progress_update(sid, "bootstrap", "start", sprintf("n_boot=%d", n_boot))
  bio_modules <- setdiff(sort(unique(mods$module)), c("grey", "gold"))
  base_assign <- setNames(mods$module, mods$taxon)
  boot_rows <- rbindlist(lapply(seq_len(n_boot), function(i) {
    if (i == 1 || i == n_boot || i %% 10 == 0) {
      log_msg(sprintf("[%s] bootstrap replicate %d/%d", sid, i, n_boot))
      progress_update(sid, "bootstrap_replicate", "ok", sprintf("%d/%d", i, n_boot))
    }
    set.seed(PARAMS$seed + i)
    sample_override <- lapply(PARAMS$stage1_cores, function(core_id) {
      sample(seq_len(nrow(expr_by_core[[core_id]])), size = nrow(expr_by_core[[core_id]]), replace = TRUE)
    })
    names(sample_override) <- PARAMS$stage1_cores
    fit <- tryCatch(
      fit_consensus(par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize, sample_override),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    asg <- setNames(as.character(fit$colors), names(fit$colors))
    asg <- asg[names(base_assign)]
    boot_mods <- setdiff(sort(unique(asg)), c("grey", "gold"))
    rbindlist(lapply(bio_modules, function(m) {
      bg <- names(base_assign)[base_assign == m]
      if (!length(boot_mods)) return(data.table(replicate = i, module = m, jaccard = NA_real_))
      js <- sapply(boot_mods, function(bm) jaccard(bg, names(asg)[asg == bm]))
      data.table(replicate = i, module = m, jaccard = max(js, na.rm = TRUE))
    }))
  }), fill = TRUE)
  fwrite(boot_rows, file.path(sdir, "bootstrap_module_stability.tsv"), sep = "\t")
  boot_sum <- boot_rows[, .(
    n = .N,
    jaccard_median = median(jaccard, na.rm = TRUE),
    jaccard_p05 = quantile(jaccard, 0.05, na.rm = TRUE),
    jaccard_p95 = quantile(jaccard, 0.95, na.rm = TRUE)
  ), by = module]
  fwrite(boot_sum, file.path(sdir, "bootstrap_module_stability_summary.tsv"), sep = "\t")
  progress_update(sid, "bootstrap", "ok")

  progress_update(sid, "core_balance", "start")
  core_sizes <- meta_filt[core %in% PARAMS$stage1_cores, .N, by = core]
  min_n <- min(core_sizes$N)
  set.seed(PARAMS$seed + 999)
  bal_idx <- lapply(PARAMS$stage1_cores, function(core_id) {
    sample(seq_len(nrow(expr_by_core[[core_id]])), size = min_n, replace = FALSE)
  })
  names(bal_idx) <- PARAMS$stage1_cores
  net_bal <- fit_consensus(par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize, bal_idx)
  bal_asg <- setNames(as.character(net_bal$colors), names(net_bal$colors))
  full_asg <- setNames(as.character(net$colors), names(net$colors))
  modules_full <- setdiff(sort(unique(full_asg)), c("grey", "gold"))
  imbalance_dt <- rbindlist(lapply(modules_full, function(m) {
    fg <- names(full_asg)[full_asg == m]
    bm <- setdiff(sort(unique(bal_asg)), c("grey", "gold"))
    if (!length(bm)) return(data.table(module = m, best_balanced_jaccard = NA_real_))
    js <- sapply(bm, function(x) jaccard(fg, names(bal_asg)[bal_asg == x]))
    data.table(module = m, best_balanced_jaccard = max(js, na.rm = TRUE))
  }))
  fwrite(imbalance_dt, file.path(sdir, "core_balance_module_jaccard.tsv"), sep = "\t")
  progress_update(sid, "core_balance", "ok")

  summary_row <- data.table(
    setting_id = sid,
    power = par$power,
    deepSplit = par$deepSplit,
    mergeCutHeight = par$mergeCutHeight,
    minModuleSize = par$minModuleSize,
    non_grey_modules = mods[module != "grey", uniqueN(module)],
    grey_pct = mods[module == "grey", .N] / nrow(mods) * 100,
    bio_pres_strong = pres[module_type == "biological" & preserved == "strong", .N],
    bio_pres_moderate = pres[module_type == "biological" & preserved == "moderate", .N],
    mean_concordance_pearson = mean(conc$pearson_r, na.rm = TRUE),
    mean_concordance_spearman = mean(conc$spearman_rho, na.rm = TRUE),
    mean_concordance_rmse = mean(conc$rmse, na.rm = TRUE),
    mean_bootstrap_jaccard = mean(boot_sum$jaccard_median, na.rm = TRUE),
    mean_balanced_jaccard = mean(imbalance_dt$best_balanced_jaccard, na.rm = TRUE)
  )
  fwrite(summary_row, file.path(sdir, "setting_summary.tsv"), sep = "\t")
  progress_update(sid, "done", "ok")
  summary_row
}

res <- rbindlist(lapply(seq_len(nrow(settings)), function(i) {
  par <- settings[i]
  tryCatch(
    evaluate_setting(par),
    error = function(e) {
      log_msg(sprintf("ERROR in %s: %s", par$setting_id, conditionMessage(e)))
      progress_update(par$setting_id, "failed", "error", conditionMessage(e))
      data.table(
        setting_id = par$setting_id,
        power = par$power, deepSplit = par$deepSplit,
        mergeCutHeight = par$mergeCutHeight, minModuleSize = par$minModuleSize,
        non_grey_modules = NA_integer_, grey_pct = NA_real_,
        bio_pres_strong = NA_integer_, bio_pres_moderate = NA_integer_,
        mean_concordance_pearson = NA_real_, mean_concordance_spearman = NA_real_,
        mean_concordance_rmse = NA_real_, mean_bootstrap_jaccard = NA_real_,
        mean_balanced_jaccard = NA_real_
      )
    }
  )
}), fill = TRUE)

fwrite(res, file.path(OUT_TABLE, "filtered_candidate_summary.tsv"), sep = "\t")
ranked <- rank_settings(res)
fwrite(ranked, file.path(OUT_TABLE, "filtered_candidate_ranked.tsv"), sep = "\t")

orig <- fread(here("networkQC", "results", "full_eval", "all_settings_ranked.tsv"))[
  setting_id %in% settings$setting_id
]
compare <- merge(
  orig[, .(
    setting_id, original_rank = rank, original_final_score = final_score,
    original_grey_pct = grey_pct, original_non_grey_modules = non_grey_modules,
    original_bootstrap = mean_bootstrap_jaccard,
    original_balance = mean_balanced_jaccard,
    original_pres_strong = bio_pres_strong,
    original_concordance = mean_concordance_pearson
  )],
  ranked[, .(
    setting_id, filtered_rank = sensitivity_rank, filtered_final_score = final_score,
    filtered_grey_pct = grey_pct, filtered_non_grey_modules = non_grey_modules,
    filtered_bootstrap = mean_bootstrap_jaccard,
    filtered_balance = mean_balanced_jaccard,
    filtered_pres_strong = bio_pres_strong,
    filtered_concordance = mean_concordance_pearson
  )],
  by = "setting_id",
  all = TRUE
)
compare[, `:=`(
  rank_change = filtered_rank - original_rank,
  delta_grey_pct = filtered_grey_pct - original_grey_pct,
  delta_non_grey_modules = filtered_non_grey_modules - original_non_grey_modules,
  delta_bootstrap = filtered_bootstrap - original_bootstrap,
  delta_balance = filtered_balance - original_balance,
  delta_pres_strong = filtered_pres_strong - original_pres_strong,
  delta_concordance = filtered_concordance - original_concordance
)]
setorder(compare, filtered_rank)
fwrite(compare, file.path(OUT_TABLE, "original_vs_filtered_candidate_comparison.tsv"), sep = "\t")

overlap_rows <- rbindlist(lapply(settings$setting_id, function(sid) {
  fpath <- file.path(OUT_SETTINGS, sid, "module_assignments.tsv")
  opath <- here("networkQC", "results", "full_eval", sid, "module_assignments.tsv")
  if (!file.exists(fpath) || !file.exists(opath)) return(NULL)
  filt <- fread(fpath)
  org <- fread(opath)
  common_taxa <- intersect(filt$taxon, org$taxon)
  filt_asg <- setNames(filt$module, filt$taxon)[common_taxa]
  org_asg <- setNames(org$module, org$taxon)[common_taxa]
  org_mods <- setdiff(sort(unique(org_asg)), c("grey", "gold"))
  filt_mods <- setdiff(sort(unique(filt_asg)), c("grey", "gold"))
  rbindlist(lapply(org_mods, function(om) {
    taxa <- names(org_asg)[org_asg == om]
    js <- sapply(filt_mods, function(fm) jaccard(taxa, names(filt_asg)[filt_asg == fm]))
    data.table(
      setting_id = sid,
      original_module = om,
      best_filtered_module = names(which.max(js)),
      best_jaccard = max(js, na.rm = TRUE),
      original_size = length(taxa),
      filtered_size = sum(filt_asg == names(which.max(js)))
    )
  }))
}), fill = TRUE)
fwrite(overlap_rows, file.path(OUT_TABLE, "original_vs_filtered_module_overlap.tsv"), sep = "\t")

overlap_summary <- overlap_rows[, .(
  mean_module_overlap = mean(best_jaccard, na.rm = TRUE),
  median_module_overlap = median(best_jaccard, na.rm = TRUE),
  min_module_overlap = min(best_jaccard, na.rm = TRUE)
), by = setting_id]
fwrite(overlap_summary, file.path(OUT_TABLE, "original_vs_filtered_module_overlap_summary.tsv"), sep = "\t")

p_rank <- ggplot(compare, aes(reorder(setting_id, filtered_final_score), filtered_final_score, fill = setting_id)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("rank %d", filtered_rank)), vjust = -0.3, size = 3) +
  coord_cartesian(ylim = c(0, max(compare$filtered_final_score, na.rm = TRUE) * 1.15)) +
  labs(title = "Candidate ranking after removing ST8 low-taxa samples", x = NULL, y = "Filtered sensitivity score") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(OUT_FIG, "filtered_candidate_ranking.png"), p_rank, width = 7, height = 4.8, dpi = 160)

heat_dt <- melt(
  compare[, .(setting_id, delta_grey_pct, delta_bootstrap, delta_balance, delta_pres_strong, delta_concordance)],
  id.vars = "setting_id",
  variable.name = "metric",
  value.name = "delta"
)
p_delta <- ggplot(heat_dt, aes(metric, setting_id, fill = delta)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", delta)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, na.value = "grey90") +
  labs(title = "Filtered minus original metric changes", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_FIG, "original_vs_filtered_metric_deltas.png"), p_delta, width = 8.5, height = 4.5, dpi = 160)

top_setting <- ranked[order(sensitivity_rank)][1, setting_id]
orig_top <- orig[order(original_rank)][1, setting_id]

sink(REPORT)
cat("# ST8 Low-Taxa Removal Network Sensitivity\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat(sprintf("- Removed samples: `%d` ST8 low-taxa samples (`detected_taxa < 250`)\n", length(remove_samples)))
cat(sprintf("- Build settings: `n_boot=%d`, `n_perm=%d`\n", n_boot, n_perm))
cat("- Progress log: `InputQC/st8_low_taxa_network_sensitivity/results/st8_low_taxa_network_sensitivity.log`\n\n")

cat("## Main finding\n\n")
if (identical(top_setting, orig_top)) {
  cat(sprintf("The top candidate did **not** change: `%s` remains the top-ranked setting after removing ST8 low-taxa samples.\n\n", top_setting))
} else {
  cat(sprintf("The top candidate changed from `%s` to `%s` after removing ST8 low-taxa samples.\n\n", orig_top, top_setting))
}

cat("## Candidate comparison\n\n")
cat("|setting|original rank|filtered rank|filtered score|grey pct|non-grey modules|bootstrap Jaccard|balanced Jaccard|strong preserved|Pearson concordance|\n")
cat("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(compare))) {
  cat(sprintf("|%s|%d|%d|%.3f|%.2f|%d|%.3f|%.3f|%d|%.3f|\n",
              compare$setting_id[i], compare$original_rank[i], compare$filtered_rank[i],
              compare$filtered_final_score[i], compare$filtered_grey_pct[i],
              compare$filtered_non_grey_modules[i], compare$filtered_bootstrap[i],
              compare$filtered_balance[i], compare$filtered_pres_strong[i],
              compare$filtered_concordance[i]))
}
cat("\n")

cat("## Module overlap with original full-data modules\n\n")
cat("|setting|mean overlap|median overlap|min overlap|\n")
cat("|---|---:|---:|---:|\n")
for (i in seq_len(nrow(overlap_summary))) {
  cat(sprintf("|%s|%.3f|%.3f|%.3f|\n",
              overlap_summary$setting_id[i], overlap_summary$mean_module_overlap[i],
              overlap_summary$median_module_overlap[i], overlap_summary$min_module_overlap[i]))
}
cat("\n")

cat("## Interpretation note\n\n")
cat("This is a build-depth sensitivity run, not an overnight final run. It is designed to detect large rank flips or instability caused by removing the ST8 low-taxa samples. If the top two candidates are close or the conclusion is borderline, rerun this same script with larger `--n_boot` and `--n_perm`.\n\n")

cat("## Outputs\n\n")
cat("- Filtered ranking: `InputQC/st8_low_taxa_network_sensitivity/results/tables/filtered_candidate_ranked.tsv`\n")
cat("- Original vs filtered comparison: `InputQC/st8_low_taxa_network_sensitivity/results/tables/original_vs_filtered_candidate_comparison.tsv`\n")
cat("- Module overlap: `InputQC/st8_low_taxa_network_sensitivity/results/tables/original_vs_filtered_module_overlap.tsv`\n")
cat("- Figures: `InputQC/st8_low_taxa_network_sensitivity/results/figures/`\n")
sink()

log_msg("Report written: ", REPORT)
log_msg("Complete")
