#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

source(here::here("balancednetwork", "config_balanced.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA_character_) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}

run_mode <- arg_val("mode", "build")
if (!run_mode %in% c("build", "final")) stop("Invalid mode.")
n_boot <- as.integer(arg_val("n_boot", ifelse(run_mode == "final", "1000", "100")))
n_perm <- as.integer(arg_val("n_perm", ifelse(run_mode == "final", "700", "200")))
top_n <- as.integer(arg_val("top_n", "5"))
force <- as.integer(arg_val("force", "0")) == 1L
setting_rank <- arg_val("setting_rank", "")
setting_id_filter <- arg_val("setting_id", "")

progress_tsv <- Sys.getenv("BAL_PROGRESS_TSV", unset = "")
if (!nzchar(progress_tsv)) {
  progress_tsv <- file.path(BAL$logs_dir, sprintf("full_eval_progress_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
if (!file.exists(progress_tsv)) {
  fwrite(
    data.table(
      timestamp = character(),
      run_mode = character(),
      setting_id = character(),
      stage = character(),
      status = character(),
      detail = character()
    ),
    progress_tsv,
    sep = "\t"
  )
}

progress_update <- function(setting_id, stage, status, detail = "") {
  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    run_mode = run_mode,
    setting_id = setting_id,
    stage = stage,
    status = status,
    detail = detail
  )
  fwrite(row, progress_tsv, sep = "\t", append = TRUE)
}

ctx <- load_balanced_context()
meta <- ctx$meta
expr_by_core <- ctx$expr_by_core
expr_validation <- ctx$expr_validation
pool_by_core_bin <- ctx$pool_by_core_bin
quotas <- ctx$quotas

sel_all <- fread(file.path(BAL$qc_tables_dir, "qc_decision_matrix.tsv"))
sel_all[, setting_id := paste0("top", rank)]
if (nzchar(setting_rank)) {
  sel <- sel_all[rank == as.integer(setting_rank)]
} else if (nzchar(setting_id_filter)) {
  sel <- sel_all[setting_id == setting_id_filter]
} else {
  sel <- sel_all[1:min(top_n, .N)]
  if (nrow(sel) > 0) {
    cutoff <- sel[nrow(sel), decision_score]
    sel <- sel_all[decision_score >= cutoff]
  }
}
if (nrow(sel) == 0) stop("No settings matched the requested selector.")
fwrite(sel, file.path(BAL$qc_full_dir, "settings_to_evaluate.tsv"), sep = "\t")
progress_update("global", "select_settings", "ok", sprintf("selected_n=%d", nrow(sel)))

r1 <- "GeoB25202_R1"
r2 <- BAL_PARAMS$validation_core

age_aligned_concordance <- function(MEs_all, me_cols) {
  meta_age <- meta[, .(label, core, age_kyr = y_bp / 1000)]
  dt <- merge(MEs_all, meta_age, by.x = "sample", by.y = "label", all.x = TRUE)
  rbindlist(lapply(me_cols, function(me) {
    d1 <- dt[core == r1, .(age_kyr, v = get(me))][order(age_kyr)]
    d2 <- dt[core == r2, .(age_kyr, v = get(me))][order(age_kyr)]
    d1 <- d1[is.finite(age_kyr) & is.finite(v)]
    d2 <- d2[is.finite(age_kyr) & is.finite(v)]
    if (nrow(d1) < 3 || nrow(d2) < 3) return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    lo <- max(min(d1$age_kyr), min(d2$age_kyr))
    hi <- min(max(d1$age_kyr), max(d2$age_kyr))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    x <- seq(lo, hi, length.out = BAL_PARAMS$age_grid_points)
    y1 <- approx(d1$age_kyr, d1$v, xout = x, rule = 2)$y
    y2 <- approx(d2$age_kyr, d2$v, xout = x, rule = 2)$y
    data.table(
      module = me,
      pearson_r = unname(cor(y1, y2, method = "pearson")),
      spearman_rho = unname(cor(y1, y2, method = "spearman")),
      rmse = sqrt(mean((y1 - y2)^2, na.rm = TRUE))
    )
  }))
}

sample_balanced_labels <- function(seed_offset, replace = TRUE) {
  set.seed(PARAMS$seed + seed_offset)
  out <- lapply(BAL_PARAMS$train_cores, function(core_id) {
    pools <- pool_by_core_bin[[core_id]]
    unlist(lapply(seq_len(nrow(quotas)), function(i) {
      bin_id <- quotas$age_bin[i]
      quota <- quotas$quota_per_core[i]
      pool <- pools[[bin_id]]
      if (is.null(pool) || quota <= 0) return(character(0))
      sample(pool, size = quota, replace = replace)
    }), use.names = FALSE)
  })
  names(out) <- BAL_PARAMS$train_cores
  out
}

results_rows <- list()
all_fail_rows <- list()

for (k in seq_len(nrow(sel))) {
  par <- sel[k]
  sid <- par$setting_id
  sdir <- file.path(BAL$qc_full_dir, sid)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  summary_path <- file.path(sdir, "setting_summary.tsv")
  progress_update(sid, "setting", "start", sprintf("power=%s deepSplit=%s mergeCutHeight=%s minModuleSize=%s", par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize))
  if (file.exists(summary_path) && !force) {
    log_msg(sprintf("[%s] exists, skipping (use --force=1 to rebuild)", sid))
    progress_update(sid, "setting", "skipped", "setting_summary.tsv exists and force=0")
    results_rows[[k]] <- fread(summary_path)
    next
  }

  log_msg(sprintf("[%s] fitting baseline balanced network", sid))
  progress_update(sid, "baseline_fit", "start")
  base_fit <- fit_consensus_from_expr(
    expr_by_core = expr_by_core,
    power = par$power,
    deepSplit = par$deepSplit,
    mergeCutHeight = par$mergeCutHeight,
    minModuleSize = par$minModuleSize
  )
  progress_update(sid, "baseline_fit", "ok")

  mods <- data.table(taxon = names(base_fit$colors), module = as.character(base_fit$colors))
  fwrite(mods, file.path(sdir, "module_assignments.tsv"), sep = "\t")
  base_assign <- setNames(mods$module, mods$taxon)
  bio_modules <- setdiff(sort(unique(mods$module)), c("grey", "gold"))

  pooled <- do.call(rbind, expr_by_core)
  MEs_train <- orderMEs(moduleEigengenes(pooled, base_fit$colors)$eigengenes)
  MEs_valid <- orderMEs(moduleEigengenes(expr_validation, base_fit$colors)$eigengenes)
  MEs_all <- rbind(
    as.data.table(MEs_train, keep.rownames = "sample"),
    as.data.table(MEs_valid, keep.rownames = "sample")
  )
  me_cols <- setdiff(names(MEs_all), "sample")
  conc <- age_aligned_concordance(MEs_all, me_cols)
  fwrite(conc, file.path(sdir, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
  progress_update(sid, "eigengene_concordance", "ok", sprintf("modules=%d", nrow(conc)))

  var_r1 <- apply(expr_by_core[[r1]], 2, var)
  var_r2 <- apply(expr_validation, 2, var)
  good <- names(which(var_r1 > 0 & var_r2 > 0))
  mp <- modulePreservation(
    multiData = list(
      GeoB_R1 = list(data = expr_by_core[[r1]][, good, drop = FALSE]),
      GeoB_R2 = list(data = expr_validation[, good, drop = FALSE])
    ),
    multiColor = list(GeoB_R1 = base_assign[good]),
    referenceNetworks = 1,
    testNetworks = 2,
    nPermutations = n_perm,
    randomSeed = PARAMS$seed,
    verbose = 0
  )
  ps <- mp$preservation$Z[[1]][[2]]
  pres <- data.table(module = rownames(ps), Zsummary = ps$Zsummary.pres)
  pres[, preserved := fcase(Zsummary > 10, "strong", Zsummary > 2, "moderate", default = "weak")]
  pres[, module_type := fcase(module %in% c("grey", "gold"), "technical", default = "biological")]
  fwrite(pres, file.path(sdir, "preservation.tsv"), sep = "\t")
  progress_update(sid, "preservation", "ok", sprintf("rows=%d n_perm=%d", nrow(pres), n_perm))

  log_msg(sprintf("[%s] running balanced bootstrap n=%d", sid, n_boot))
  progress_update(sid, "bootstrap", "start", sprintf("n_boot=%d", n_boot))
  boot_rows <- vector("list", n_boot)
  fail_rows <- list()
  for (i in seq_len(n_boot)) {
    if (i %% 25 == 0 || i == 1 || i == n_boot) {
      log_msg(sprintf("[%s] bootstrap %d/%d", sid, i, n_boot))
      progress_update(sid, "bootstrap_progress", "ok", sprintf("%d/%d", i, n_boot))
    }
    labels <- sample_balanced_labels(seed_offset = (k * 100000 + i), replace = TRUE)
    expr_boot <- lapply(BAL_PARAMS$train_cores, function(core_id) {
      ctx$vst[labels[[core_id]], names(base_assign), drop = FALSE]
    })
    names(expr_boot) <- BAL_PARAMS$train_cores
    fit <- tryCatch(
      fit_consensus_from_expr(
        expr_by_core = expr_boot,
        power = par$power,
        deepSplit = par$deepSplit,
        mergeCutHeight = par$mergeCutHeight,
        minModuleSize = par$minModuleSize
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      fail_rows[[length(fail_rows) + 1L]] <- data.table(replicate = i, error = conditionMessage(fit))
      progress_update(sid, "bootstrap_replicate", "failed", sprintf("replicate=%d error=%s", i, conditionMessage(fit)))
      next
    }
    asg <- setNames(as.character(fit$colors), names(fit$colors))
    asg <- asg[names(base_assign)]
    boot_mods <- setdiff(sort(unique(asg)), c("grey", "gold"))
    boot_rows[[i]] <- rbindlist(lapply(bio_modules, function(m) {
      bg <- names(base_assign)[base_assign == m]
      if (!length(boot_mods)) return(data.table(replicate = i, module = m, jaccard = NA_real_))
      js <- sapply(boot_mods, function(bm) jaccard(bg, names(asg)[asg == bm]))
      data.table(replicate = i, module = m, jaccard = max(js, na.rm = TRUE))
    }))
  }
  boot_dt <- rbindlist(boot_rows, fill = TRUE)
  fail_dt <- if (length(fail_rows)) rbindlist(fail_rows) else data.table(replicate = integer(), error = character())
  if (nrow(boot_dt) == 0) stop(sprintf("[%s] no successful bootstrap results", sid))
  boot_sum <- boot_dt[, .(
    n = .N,
    jaccard_median = median(jaccard, na.rm = TRUE),
    jaccard_p05 = quantile(jaccard, 0.05, na.rm = TRUE),
    jaccard_p95 = quantile(jaccard, 0.95, na.rm = TRUE)
  ), by = module]
  fwrite(boot_dt, file.path(sdir, "bootstrap_module_stability.tsv"), sep = "\t")
  fwrite(boot_sum, file.path(sdir, "bootstrap_module_stability_summary.tsv"), sep = "\t")
  fwrite(fail_dt, file.path(sdir, "bootstrap_failures.tsv"), sep = "\t")
  progress_update(sid, "bootstrap", "ok", sprintf("successful=%d failed=%d", uniqueN(boot_dt$replicate), nrow(fail_dt)))
  if (nrow(fail_dt) > 0) all_fail_rows[[length(all_fail_rows) + 1L]] <- cbind(data.table(setting_id = sid), fail_dt)

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
    n_boot_requested = n_boot,
    n_boot_successful = uniqueN(boot_dt$replicate),
    n_boot_failed = nrow(fail_dt)
  )
  fwrite(summary_row, summary_path, sep = "\t")
  results_rows[[k]] <- summary_row
  progress_update(sid, "setting", "ok", sprintf("final_score_inputs boot=%.4f concordance=%.4f", summary_row$mean_bootstrap_jaccard, summary_row$mean_concordance_pearson))
}

res <- rbindlist(results_rows, fill = TRUE)
score_input <- copy(res)
if (nzchar(setting_rank) || nzchar(setting_id_filter)) {
  summary_dirs <- list.dirs(BAL$qc_full_dir, recursive = FALSE, full.names = TRUE)
  summary_dirs <- summary_dirs[grepl("/top[0-9]+$", summary_dirs)]
  summary_rows <- rbindlist(lapply(summary_dirs, function(d) {
    f <- file.path(d, "setting_summary.tsv")
    if (file.exists(f)) fread(f) else NULL
  }), fill = TRUE)
  if (nrow(summary_rows) > 0) score_input <- summary_rows
}
score <- copy(score_input)
score[, s_grey := norm01(grey_pct, FALSE)]
score[, s_boot := norm01(mean_bootstrap_jaccard, TRUE)]
score[, s_pres := norm01(bio_pres_strong + 0.5 * bio_pres_moderate, TRUE)]
score[, s_conc := norm01(mean_concordance_pearson, TRUE)]
score[, final_score := 0.30 * s_grey + 0.35 * s_boot + 0.20 * s_pres + 0.15 * s_conc]
setorder(score, -final_score, grey_pct)
score[, rank := .I]

fwrite(score, file.path(BAL$qc_full_dir, "all_settings_ranked.tsv"), sep = "\t")
if (length(all_fail_rows)) fwrite(rbindlist(all_fail_rows), file.path(BAL$qc_full_dir, "all_bootstrap_failures.tsv"), sep = "\t")
progress_update("global", "ranking", "ok", sprintf("ranked_n=%d", nrow(score)))

log_msg("balanced full evaluation complete")
progress_update("global", "run", "ok", "balanced full evaluation complete")
