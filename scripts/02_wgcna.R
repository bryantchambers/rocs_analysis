#!/usr/bin/env Rscript
# 02_wgcna.R — Consensus WGCNA prokaryote co-expression modules
#
# Consensus across ST8, ST13, GeoB25202_R1.
# GeoB25202_R2 held out for preservation testing.
# Soft power selected by signed R² (negative slope required).
#
# Outputs: results/stage1/wgcna/
#   module_assignments.tsv    taxon → module colour
#   module_eigengenes.tsv     sample × module eigengenes (all cores)
#   soft_power.tsv            scale-free fit table
#   preservation.tsv          Zsummary per module (R1 → R2)
#   consensus_wgcna.rds       full blockwiseConsensusModules object (for inspection)

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
preservation_perms <- if (run_mode == "final") {
  PARAMS$wgcna_preservation_permutations_final
} else {
  PARAMS$wgcna_preservation_permutations_build
}
log_msg(sprintf("Run mode: %s (modulePreservation permutations=%d)", run_mode, preservation_perms))

WGCNA_OUT <- file.path(RESULTS$stage1, "wgcna")
dir.create(WGCNA_OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load ───────────────────────────────────────────────────────────────────

log_msg("Loading VST matrix and metadata...")
vst  <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))

# Split by core
expr_by_core <- lapply(PARAMS$all_cores, function(core_id) {
  samps <- meta[core == core_id, label]
  samps <- intersect(samps, rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- PARAMS$all_cores

for (c in PARAMS$all_cores)
  log_msg(sprintf("  %s: %d samples", c, nrow(expr_by_core[[c]])))

# ── 2. Soft-power selection (signed R²: negative slope required) ──────────────

log_msg("Selecting soft-thresholding power (signed R²)...")

powers <- c(1:10, seq(12, 20, 2))

sft_list <- lapply(PARAMS$stage1_cores, function(core_id) {
  sft <- pickSoftThreshold(expr_by_core[[core_id]],
                           powerVector = powers,
                           networkType = "signed",
                           verbose = 0)
  dt <- as.data.table(sft$fitIndices)
  dt[, core := core_id]
  dt
})
sft_dt <- rbindlist(sft_list)
fwrite(sft_dt, file.path(WGCNA_OUT, "soft_power.tsv"), sep = "\t")

# signedR2: only credit negative slope (proper scale-free topology)
sft_dt[, signedR2 := ifelse(slope < 0, SFT.R.sq, 0)]

# Per-core passing powers (signedR2 ≥ 0.80)
passing <- sft_dt[signedR2 >= 0.80, .(min_pass = min(Power)), by = core]

if (nrow(passing) == length(PARAMS$stage1_cores)) {
  soft_power <- max(passing$min_pass)
  log_msg(sprintf("  All cores pass R²≥0.80 — using power %d (max across cores)", soft_power))
} else {
  # Fallback: best signed R² with negative slope and mean.k. > 5.
  # mean.k. > 5 ensures the network is not too sparse (< 5 expected connections
  # per node would make module detection unreliable).
  best <- sft_dt[slope < 0 & mean.k. > 5, .(best_r2 = max(signedR2)), by = Power]
  soft_power <- if (nrow(best) > 0) best[which.max(best_r2), Power] else 12L
  log_msg(sprintf("  Fallback: using power %d (best signed R², mean.k.>5)", soft_power))
}

log_msg(sprintf("  Final soft power: %d", soft_power))

# ── 3. Consensus WGCNA ────────────────────────────────────────────────────────

log_msg("Building consensus WGCNA (ST8, ST13, GeoB R1)...")

# multiExpr: list of list(data=matrix) — one entry per training core
multiExpr <- lapply(PARAMS$stage1_cores, function(core_id) {
  list(data = expr_by_core[[core_id]])
})
names(multiExpr) <- PARAMS$stage1_cores

# corType = "pearson": CLR-transformed data is approximately Gaussian with zero
# mean per taxon, so Pearson correlations are appropriate and match the original
# analysis. Bicor would be more robust to outliers but is not used here.
#
# maxBlockSize = 5000: taxa are processed in blocks of ≤ 5000 to bound memory
# usage. With ~1300 taxa all fit in a single block, so this has no effect on
# results here. Increase if taxon count grows beyond 5000.
net <- blockwiseConsensusModules(
  multiExpr       = multiExpr,
  power           = soft_power,
  networkType     = "signed",
  corType         = "pearson",
  maxBlockSize    = 5000,
  minModuleSize   = PARAMS$wgcna_min_module_size,
  deepSplit       = PARAMS$wgcna_deep_split,
  mergeCutHeight  = PARAMS$wgcna_merge_cut_height,
  numericLabels   = FALSE,
  saveTOMs        = FALSE,
  verbose         = 2
)

module_colors <- net$colors
n_mod  <- length(unique(module_colors)) - 1
n_grey <- sum(module_colors == "grey")

log_msg(sprintf("  Modules: %d non-grey, %d grey taxa", n_mod, n_grey))
for (m in sort(names(table(module_colors))))
  log_msg(sprintf("    %s: %d", m, sum(module_colors == m)))

# ── 4. Module eigengenes (pooled training + validation) ───────────────────────

log_msg("Computing module eigengenes...")

# Pool all samples for eigengene computation
pooled <- do.call(rbind, lapply(PARAMS$stage1_cores, function(c) expr_by_core[[c]]))
MEs_train <- moduleEigengenes(pooled, module_colors)$eigengenes
MEs_train <- orderMEs(MEs_train)

# Validation core eigengenes
MEs_valid <- moduleEigengenes(expr_by_core[[PARAMS$validation_core]],
                               module_colors)$eigengenes
MEs_valid <- orderMEs(MEs_valid)

MEs_all <- rbind(
  as.data.table(MEs_train, keep.rownames = "sample"),
  as.data.table(MEs_valid, keep.rownames = "sample")
)

# R1/R2 concordance
me_cols  <- setdiff(names(MEs_all), "sample")
r1_samps <- rownames(expr_by_core[["GeoB25202_R1"]])
r2_samps <- rownames(expr_by_core[[PARAMS$validation_core]])

# R1 and R2 have different sample counts (26 vs 25), so a direct vector
# correlation is not meaningful. We log per-module means as a sanity check;
# for a formal concordance test, interpolate both cores to a common age grid.
log_msg("  R1/R2 eigengene mean (quick sanity check)")
for (me in me_cols) {
  m1 <- mean(MEs_all[sample %in% r1_samps][[me]], na.rm = TRUE)
  m2 <- mean(MEs_all[sample %in% r2_samps][[me]], na.rm = TRUE)
  log_msg(sprintf("    %-12s R1_mean=%.3f  R2_mean=%.3f", me, m1, m2))
}

# ── 5. Module preservation (R1 → R2) ─────────────────────────────────────────

log_msg("Module preservation (GeoB R1 → R2)...")

var_r1 <- apply(expr_by_core[["GeoB25202_R1"]], 2, var)
var_r2 <- apply(expr_by_core[[PARAMS$validation_core]], 2, var)
good   <- names(which(var_r1 > 0 & var_r2 > 0))

mp <- modulePreservation(
  multiData = list(
    GeoB_R1 = list(data = expr_by_core[["GeoB25202_R1"]][, good]),
    GeoB_R2 = list(data = expr_by_core[[PARAMS$validation_core]][, good])
  ),
  multiColor       = list(GeoB_R1 = module_colors[good]),
  referenceNetworks = 1,
  testNetworks     = 2,
  nPermutations    = preservation_perms,
  randomSeed       = PARAMS$seed,
  verbose          = 0
)

pres_stats <- mp$preservation$Z[[1]][[2]]
pres_dt <- data.table(
  module        = rownames(pres_stats),
  Zsummary      = pres_stats$Zsummary.pres,
  Zdensity      = pres_stats$Zdensity.pres,
  Zconnectivity = pres_stats$Zconnectivity.pres
)
pres_dt[, preserved := fcase(
  Zsummary > 10, "strong",
  Zsummary >  2, "moderate",
  default = "weak"
)]
pres_dt[, module_type := fcase(
  module %in% c("grey", "gold"), "technical",
  default = "biological"
)]
log_msg("  Preservation:")
for (p in c("strong", "moderate", "weak"))
  log_msg(sprintf("    %s: %d modules", p, sum(pres_dt$preserved == p, na.rm = TRUE)))

pres_bio <- pres_dt[module_type == "biological"]
log_msg("  Biological module preservation only:")
for (p in c("strong", "moderate", "weak"))
  log_msg(sprintf("    %s: %d modules", p, sum(pres_bio$preserved == p, na.rm = TRUE)))

# ── 5b. Age-aligned eigengene concordance (R1 ↔ R2) ──────────────────────────

log_msg("Age-aligned R1/R2 eigengene concordance on common age grid...")
r1_dt <- merge(
  data.table(sample = r1_samps),
  meta[, .(label, y_bp)],
  by.x = "sample",
  by.y = "label",
  all.x = TRUE
)
r2_dt <- merge(
  data.table(sample = r2_samps),
  meta[, .(label, y_bp)],
  by.x = "sample",
  by.y = "label",
  all.x = TRUE
)
r1_dt[, age_kyr := y_bp / 1000]
r2_dt[, age_kyr := y_bp / 1000]

me_train <- MEs_all[sample %in% r1_samps]
me_valid <- MEs_all[sample %in% r2_samps]

concordance <- rbindlist(lapply(me_cols, function(me) {
  d1 <- merge(r1_dt[, .(sample, age_kyr)], me_train[, .(sample, value = get(me))], by = "sample")
  d2 <- merge(r2_dt[, .(sample, age_kyr)], me_valid[, .(sample, value = get(me))], by = "sample")
  d1 <- d1[is.finite(age_kyr) & is.finite(value)][order(age_kyr)]
  d2 <- d2[is.finite(age_kyr) & is.finite(value)][order(age_kyr)]
  if (nrow(d1) < 3 || nrow(d2) < 3) {
    return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
  }
  age_min <- max(min(d1$age_kyr), min(d2$age_kyr))
  age_max <- min(max(d1$age_kyr), max(d2$age_kyr))
  if (!is.finite(age_min) || !is.finite(age_max) || age_max <= age_min) {
    return(data.table(module = me, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_))
  }
  xout <- seq(age_min, age_max, length.out = PARAMS$wgcna_stability_age_grid_points)
  y1 <- approx(d1$age_kyr, d1$value, xout = xout, rule = 2)$y
  y2 <- approx(d2$age_kyr, d2$value, xout = xout, rule = 2)$y
  data.table(
    module = me,
    pearson_r = unname(suppressWarnings(cor(y1, y2, method = "pearson"))),
    spearman_rho = unname(suppressWarnings(cor(y1, y2, method = "spearman"))),
    rmse = sqrt(mean((y1 - y2)^2, na.rm = TRUE))
  )
}), fill = TRUE)
setnames(concordance,
         old = names(concordance),
         new = sub("\\.V1$", "", names(concordance)))

# ── 6. Save ───────────────────────────────────────────────────────────────────

log_msg("Saving outputs...")

fwrite(data.table(taxon = names(module_colors), module = module_colors),
       file.path(WGCNA_OUT, "module_assignments.tsv"), sep = "\t")
fwrite(MEs_all,  file.path(WGCNA_OUT, "module_eigengenes.tsv"),  sep = "\t")
fwrite(pres_dt,  file.path(WGCNA_OUT, "preservation.tsv"),       sep = "\t")
fwrite(pres_bio, file.path(WGCNA_OUT, "preservation_biological.tsv"), sep = "\t")
fwrite(concordance, file.path(WGCNA_OUT, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
saveRDS(net,     file.path(WGCNA_OUT, "consensus_wgcna.rds"))

log_msg("Done. Outputs in ", WGCNA_OUT)
