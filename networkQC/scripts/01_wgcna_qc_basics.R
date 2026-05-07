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

OUT_TABLE <- here("networkQC", "results", "tables")
OUT_FIG <- here("networkQC", "results", "figures")
dir.create(OUT_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
mods <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
me_all <- fread(file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv"))

expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- intersect(meta[core == core_id, label], rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

# 1) soft threshold curves per core
powers <- c(1:10, seq(12, 20, 2))
soft_list <- lapply(PARAMS$stage1_cores, function(core_id) {
  sft <- pickSoftThreshold(expr_by_core[[core_id]], powerVector = powers, networkType = "signed", verbose = 0)
  dt <- as.data.table(sft$fitIndices)
  dt[, core := core_id]
  dt[, signedR2 := ifelse(slope < 0, SFT.R.sq, 0)]
  dt
})
soft_dt <- rbindlist(soft_list)
fwrite(soft_dt, file.path(OUT_TABLE, "qc_soft_threshold_curves.tsv"), sep = "\t")

p1 <- ggplot(soft_dt, aes(Power, signedR2, color = core)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.8, linetype = 2) +
  theme_bw(base_size = 11) + labs(title = "Signed R2 vs Power", y = "signedR2")
ggsave(file.path(OUT_FIG, "qc_soft_threshold_signedR2.png"), p1, width = 8, height = 4, dpi = 150)

p2 <- ggplot(soft_dt, aes(Power, mean.k., color = core)) +
  geom_line() + geom_point() +
  theme_bw(base_size = 11) + labs(title = "Mean Connectivity vs Power", y = "mean.k.")
ggsave(file.path(OUT_FIG, "qc_soft_threshold_connectivity.png"), p2, width = 8, height = 4, dpi = 150)

# 2) sample clustering by core
train_samples <- unlist(lapply(PARAMS$stage1_cores, function(c) rownames(expr_by_core[[c]])))
train_mat <- vst[train_samples, , drop = FALSE]
sample_tree <- hclust(dist(train_mat), method = "average")
png(file.path(OUT_FIG, "qc_sample_dendrogram_train.png"), width = 1400, height = 800, res = 140)
plot(sample_tree, main = "Sample Dendrogram (Training Cores)", xlab = "", sub = "")
dev.off()

# 3) module size and grey burden
mod_counts <- mods[, .N, by = module][order(-N)]
mod_counts[, pct := N / sum(N) * 100]
fwrite(mod_counts, file.path(OUT_TABLE, "qc_module_size_distribution.tsv"), sep = "\t")

# 4) age-aligned eigengene concordance (train core pairs + R1/R2)
meta_age <- meta[, .(label, core, age_kyr = y_bp / 1000)]
dt_me <- merge(me_all, meta_age, by.x = "sample", by.y = "label", all.x = TRUE)
me_cols <- grep("^ME", names(dt_me), value = TRUE)

pairwise_conc <- rbindlist(lapply(combn(PARAMS$all_cores, 2, simplify = FALSE), function(pair) {
  c1 <- pair[[1]]; c2 <- pair[[2]]
  rbindlist(lapply(me_cols, function(me) {
    d1 <- dt_me[core == c1, .(age_kyr, v = get(me))][order(age_kyr)]
    d2 <- dt_me[core == c2, .(age_kyr, v = get(me))][order(age_kyr)]
    d1 <- d1[is.finite(age_kyr) & is.finite(v)]
    d2 <- d2[is.finite(age_kyr) & is.finite(v)]
    if (nrow(d1) < 3 || nrow(d2) < 3) {
      return(data.table(core_a = c1, core_b = c2, module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    }
    lo <- max(min(d1$age_kyr), min(d2$age_kyr))
    hi <- min(max(d1$age_kyr), max(d2$age_kyr))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      return(data.table(core_a = c1, core_b = c2, module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
    }
    x <- seq(lo, hi, length.out = 100)
    y1 <- approx(d1$age_kyr, d1$v, xout = x, rule = 2)$y
    y2 <- approx(d2$age_kyr, d2$v, xout = x, rule = 2)$y
    data.table(
      core_a = c1, core_b = c2, module = me,
      pearson_r = unname(cor(y1, y2, method = "pearson")),
      spearman_rho = unname(cor(y1, y2, method = "spearman")),
      rmse = sqrt(mean((y1 - y2)^2))
    )
  }))
}))
fwrite(pairwise_conc, file.path(OUT_TABLE, "qc_pairwise_core_eigengene_concordance.tsv"), sep = "\t")

# 5) whole-network quick visualization (thresholded TOM)
train_combined <- do.call(rbind, lapply(PARAMS$stage1_cores, function(c) expr_by_core[[c]]))
adj <- adjacency(train_combined, power = 20, type = "signed")
tom <- TOMsimilarityFromExpr(train_combined, networkType = "signed", TOMType = "signed", power = 20, verbose = 0)
diag(tom) <- 0

qs <- quantile(tom[upper.tri(tom)], probs = 0.995, na.rm = TRUE)
idx <- which(tom >= qs, arr.ind = TRUE)
idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
edges <- data.table(
  from = colnames(train_combined)[idx[, 1]],
  to = colnames(train_combined)[idx[, 2]],
  w = tom[idx]
)
fwrite(edges, file.path(OUT_TABLE, "qc_network_edges_top0.5pct.tsv"), sep = "\t")

summary_dt <- data.table(
  metric = c("n_samples_train", "n_taxa", "n_modules_non_grey", "grey_pct", "top_edge_threshold"),
  value = c(
    nrow(train_combined),
    ncol(train_combined),
    mods[module != "grey", uniqueN(module)],
    mods[module == "grey", .N] / nrow(mods) * 100,
    qs
  )
)
fwrite(summary_dt, file.path(OUT_TABLE, "qc_basics_summary.tsv"), sep = "\t")
log_msg("WGCNA basics QC complete")

