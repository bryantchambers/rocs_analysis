#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
  library(igraph)
})

source(here("config.R"))
set.seed(PARAMS$seed)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

OUT_BASE <- here("networkQC", "results")
OUT_TABLE <- file.path(OUT_BASE, "tables")
OUT_FIG <- file.path(OUT_BASE, "figures")
OUT_FULL <- file.path(OUT_BASE, "full_eval")
OUT_SWEEP <- file.path(OUT_BASE, "sweeps", "original")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FULL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SWEEP, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
log_file <- file.path(OUT_FULL, "progress.log")
progress_tsv <- file.path(OUT_FULL, "progress.tsv")
log_append <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", file = log_file, append = TRUE)
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
  fwrite(row, progress_tsv, sep = "\t", append = file.exists(progress_tsv), col.names = !file.exists(progress_tsv))
}

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA_character_) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}

n_boot <- as.integer(arg_val("n_boot", as.character(PARAMS$wgcna_stability_bootstrap_final)))
n_perm <- as.integer(arg_val("n_perm", as.character(PARAMS$wgcna_preservation_permutations_final)))
top_n <- as.integer(arg_val("top_n", "5"))
edge_quantile <- as.numeric(arg_val("edge_q", "0.995"))

log_msg(sprintf("full_eval settings: n_boot=%d, n_perm=%d, top_n=%d", n_boot, n_perm, top_n))
if (file.exists(log_file)) file.remove(log_file)
if (file.exists(progress_tsv)) file.remove(progress_tsv)
log_append(sprintf("full_eval settings: n_boot=%d, n_perm=%d, top_n=%d", n_boot, n_perm, top_n))

# Archive sweep artifacts in dedicated folder
for (f in c("qc_parameter_sweep_summary.tsv", "qc_parameter_sweep_recommended.tsv", "qc_decision_matrix.tsv", "qc_decision_top3.tsv")) {
  src <- file.path(OUT_TABLE, f)
  if (file.exists(src)) file.copy(src, file.path(OUT_SWEEP, f), overwrite = TRUE)
}

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
leiden_mods <- fread(file.path(OUT_TABLE, "leiden_module_assignments.tsv"))
dec <- fread(file.path(OUT_TABLE, "qc_decision_matrix.tsv"))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

train_expr <- do.call(rbind, lapply(PARAMS$stage1_cores, function(c) expr_by_core[[c]]))
all_expr <- do.call(rbind, expr_by_core)

selected <- unique(dec[, .(power, deepSplit, mergeCutHeight, minModuleSize)])[1:min(top_n, .N)]
selected[, setting_id := paste0("opt", .I)]
baseline <- data.table(
  power = 20L,
  deepSplit = PARAMS$wgcna_deep_split,
  mergeCutHeight = PARAMS$wgcna_merge_cut_height,
  minModuleSize = PARAMS$wgcna_min_module_size,
  setting_id = "baseline"
)
settings <- rbindlist(list(baseline, selected), fill = TRUE, use.names = TRUE)

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

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

age_aligned_concordance <- function(MEs_all, me_cols) {
  meta_age <- meta[, .(label, core, age_kyr = y_bp / 1000)]
  dt <- merge(MEs_all, meta_age, by.x = "sample", by.y = "label", all.x = TRUE)
  r1 <- PARAMS$stage1_cores[[length(PARAMS$stage1_cores)]]
  r2 <- PARAMS$validation_core
  rbindlist(lapply(me_cols, function(me) {
    d1 <- dt[core == r1, .(age_kyr, v = get(me))][order(age_kyr)]
    d2 <- dt[core == r2, .(age_kyr, v = get(me))][order(age_kyr)]
    d1 <- d1[is.finite(age_kyr) & is.finite(v)]
    d2 <- d2[is.finite(age_kyr) & is.finite(v)]
    if (nrow(d1) < 3 || nrow(d2) < 3) return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    lo <- max(min(d1$age_kyr), min(d2$age_kyr))
    hi <- min(max(d1$age_kyr), max(d2$age_kyr))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
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

plot_network <- function(module_assign, power, setting_label, out_png) {
  ad <- adjacency(train_expr, power = power, type = "signed")
  tom <- TOMsimilarity(ad, TOMType = "signed")
  diag(tom) <- 0
  thr <- quantile(tom[upper.tri(tom)], probs = edge_quantile, na.rm = TRUE)
  idx <- which(tom >= thr, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  edges <- data.table(from = colnames(train_expr)[idx[, 1]], to = colnames(train_expr)[idx[, 2]], w = tom[idx])
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = colnames(train_expr)))
  mods <- setNames(module_assign$module, module_assign$taxon)
  V(g)$module <- mods[V(g)$name]
  lev <- sort(unique(V(g)$module))
  pal <- setNames(rainbow(length(lev)), lev)
  colv <- pal[V(g)$module]
  png(out_png, width = 1600, height = 1200, res = 150)
  plot(g, vertex.size = 2, vertex.label = NA, edge.width = 0.2, edge.color = rgb(0, 0, 0, 0.08),
       vertex.color = colv, main = paste0("Network Graph - ", setting_label))
  legend("topleft", legend = lev, col = pal, pch = 16, pt.cex = 1, bty = "n", cex = 0.8)
  dev.off()
}

evaluate_setting <- function(par) {
  sid <- par$setting_id
  sdir <- file.path(OUT_FULL, sid)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  log_append(sprintf("Evaluating %s", sid))
  progress_update(sid, "start", "ok")

  progress_update(sid, "fit_consensus", "start")
  net <- fit_consensus(par$power, par$deepSplit, par$mergeCutHeight, par$minModuleSize)
  progress_update(sid, "fit_consensus", "ok")
  mods <- data.table(taxon = names(net$colors), module = as.character(net$colors))
  fwrite(mods, file.path(sdir, "module_assignments.tsv"), sep = "\t")
  progress_update(sid, "write_module_assignments", "ok")

  # Eigengenes for all cores
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

  # Preservation
  progress_update(sid, "preservation", "start", sprintf("n_perm=%d", n_perm))
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
  fwrite(pres, file.path(sdir, "preservation.tsv"), sep = "\t")
  fwrite(pres[module_type == "biological"], file.path(sdir, "preservation_biological.tsv"), sep = "\t")
  progress_update(sid, "preservation", "ok")

  # Age-aligned concordance
  progress_update(sid, "concordance", "start")
  me_cols <- setdiff(names(MEs_all), "sample")
  conc <- age_aligned_concordance(MEs_all, me_cols)
  fwrite(conc, file.path(sdir, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
  progress_update(sid, "concordance", "ok")

  # Bootstrap stability
  progress_update(sid, "bootstrap", "start", sprintf("n_boot=%d", n_boot))
  bio_modules <- setdiff(sort(unique(mods$module)), c("grey", "gold"))
  base_assign <- setNames(mods$module, mods$taxon)
  boot_rows <- rbindlist(lapply(seq_len(n_boot), function(i) {
    if (i %% 10 == 0 || i == 1 || i == n_boot) {
      log_append(sprintf("[%s] bootstrap replicate %d/%d", sid, i, n_boot))
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

  # Core-imbalance weighting analysis
  progress_update(sid, "core_balance", "start")
  core_sizes <- meta[core %in% PARAMS$stage1_cores, .N, by = core]
  min_n <- min(core_sizes$N)
  set.seed(PARAMS$seed + 999)
  bal_idx <- lapply(PARAMS$stage1_cores, function(core_id) sample(seq_len(nrow(expr_by_core[[core_id]])), size = min_n, replace = FALSE))
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
  imbalance_summary <- data.table(
    metric = c("train_min_core_n", "mean_best_balanced_jaccard"),
    value = c(min_n, mean(imbalance_dt$best_balanced_jaccard, na.rm = TRUE))
  )
  fwrite(imbalance_dt, file.path(sdir, "core_balance_module_jaccard.tsv"), sep = "\t")
  fwrite(imbalance_summary, file.path(sdir, "core_balance_summary.tsv"), sep = "\t")
  progress_update(sid, "core_balance", "ok")

  # Summary row
  progress_update(sid, "summary", "start")
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
  progress_update(sid, "summary", "ok")

  # Graph plot
  progress_update(sid, "plot_network", "start")
  plot_network(mods, par$power, sid, file.path(OUT_FIG, paste0("full_graph_", sid, ".png")))
  progress_update(sid, "plot_network", "ok")

  progress_update(sid, "done", "ok")
  summary_row
}

res_list <- vector("list", nrow(settings))
for (i in seq_len(nrow(settings))) {
  par <- settings[i]
  sid <- par$setting_id
  res_list[[i]] <- tryCatch(
    evaluate_setting(par),
    error = function(e) {
      msg <- paste("ERROR:", conditionMessage(e))
      log_append(sprintf("Setting %s failed: %s", sid, msg))
      progress_update(sid, "failed", "error", msg)
      data.table(
        setting_id = sid,
        power = par$power, deepSplit = par$deepSplit, mergeCutHeight = par$mergeCutHeight, minModuleSize = par$minModuleSize,
        non_grey_modules = NA_integer_, grey_pct = NA_real_,
        bio_pres_strong = NA_integer_, bio_pres_moderate = NA_integer_,
        mean_concordance_pearson = NA_real_, mean_concordance_spearman = NA_real_, mean_concordance_rmse = NA_real_,
        mean_bootstrap_jaccard = NA_real_, mean_balanced_jaccard = NA_real_
      )
    }
  )
}
res <- rbindlist(res_list, fill = TRUE)
fwrite(res, file.path(OUT_FULL, "all_settings_summary.tsv"), sep = "\t")

# Simple ranking
score <- copy(res[is.finite(grey_pct) & is.finite(mean_bootstrap_jaccard)])
norm <- function(x, hb = TRUE) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0.5, length(x)))
  z <- (x - r[1]) / (r[2] - r[1]); if (!hb) z <- 1 - z; z
}
score[, s_grey := norm(grey_pct, FALSE)]
score[, s_boot := norm(mean_bootstrap_jaccard, TRUE)]
score[, s_pres := norm(bio_pres_strong + 0.5 * bio_pres_moderate, TRUE)]
score[, s_conc := norm(mean_concordance_pearson, TRUE)]
score[, s_balance := norm(mean_balanced_jaccard, TRUE)]
score[, final_score := 0.30 * s_boot + 0.25 * s_pres + 0.20 * s_grey + 0.15 * s_conc + 0.10 * s_balance]
if (nrow(score) > 0) {
  setorder(score, -final_score)
  score[, rank := .I]
  fwrite(score, file.path(OUT_FULL, "all_settings_ranked.tsv"), sep = "\t")
} else {
  fwrite(data.table(), file.path(OUT_FULL, "all_settings_ranked.tsv"), sep = "\t")
}

# Leiden graph plot for comparison
if (file.exists(file.path(OUT_TABLE, "leiden_module_assignments.tsv"))) {
  progress_update("leiden", "plot_network", "start")
  plot_network(leiden_mods, power = unique(settings$power)[1], setting_label = "leiden", out_png = file.path(OUT_FIG, "full_graph_leiden.png"))
  progress_update("leiden", "plot_network", "ok")
}

log_append("Full eval complete")
