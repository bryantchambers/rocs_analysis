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
n_perm <- as.integer(arg_val("n_perm", ifelse(run_mode == "final", "700", "200")))
force <- as.integer(arg_val("force", "0")) == 1L

OUT_BASE <- here("balancednetwork", "results")
OUT_WGCNA <- file.path(OUT_BASE, "wgcna")
OUT_TAB <- file.path(OUT_BASE, "tables")
dir.create(OUT_WGCNA, recursive = TRUE, showWarnings = FALSE)

assign_file <- file.path(OUT_WGCNA, "module_assignments.tsv")
if (file.exists(assign_file) && !force) {
  message("[balancednetwork] balanced WGCNA exists; skipping (use --force=1 to rebuild).")
  quit(save = "no")
}

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))
sel <- fread(file.path(OUT_TAB, "balanced_baseline_samples.tsv"))

train_cores <- PARAMS$stage1_cores
valid_core <- PARAMS$validation_core

if (valid_core %in% sel$core) stop("Validation core found in training sample manifest.")

expr_by_core <- lapply(train_cores, function(core_id) {
  samps <- sel[core == core_id, label]
  samps <- intersect(samps, rownames(vst))
  vst[samps, , drop = FALSE]
})
names(expr_by_core) <- train_cores

core_n <- vapply(expr_by_core, nrow, integer(1))
if (length(unique(core_n)) != 1L) stop("Training cores are not balanced in selected sample counts.")

multiExpr <- lapply(train_cores, function(core_id) list(data = expr_by_core[[core_id]]))
names(multiExpr) <- train_cores

net <- blockwiseConsensusModules(
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
  verbose = 2
)

module_colors <- setNames(net$colors, names(net$colors))

valid_samples <- intersect(meta[core == valid_core, label], rownames(vst))
expr_valid <- vst[valid_samples, , drop = FALSE]

pooled_train <- do.call(rbind, expr_by_core)
MEs_train <- orderMEs(moduleEigengenes(pooled_train, module_colors)$eigengenes)
MEs_valid <- orderMEs(moduleEigengenes(expr_valid, module_colors)$eigengenes)
MEs_all <- rbind(
  as.data.table(MEs_train, keep.rownames = "sample"),
  as.data.table(MEs_valid, keep.rownames = "sample")
)

r1 <- "GeoB25202_R1"
r2 <- valid_core
var_r1 <- apply(expr_by_core[[r1]], 2, var)
var_r2 <- apply(expr_valid, 2, var)
good <- names(which(var_r1 > 0 & var_r2 > 0))

mp <- modulePreservation(
  multiData = list(
    GeoB_R1 = list(data = expr_by_core[[r1]][, good, drop = FALSE]),
    GeoB_R2 = list(data = expr_valid[, good, drop = FALSE])
  ),
  multiColor = list(GeoB_R1 = module_colors[good]),
  referenceNetworks = 1,
  testNetworks = 2,
  nPermutations = n_perm,
  randomSeed = PARAMS$seed,
  verbose = 0
)

pres_stats <- mp$preservation$Z[[1]][[2]]
pres_dt <- data.table(
  module = rownames(pres_stats),
  Zsummary = pres_stats$Zsummary.pres,
  Zdensity = pres_stats$Zdensity.pres,
  Zconnectivity = pres_stats$Zconnectivity.pres
)
pres_dt[, preserved := fcase(
  Zsummary > 10, "strong",
  Zsummary > 2, "moderate",
  default = "weak"
)]
pres_dt[, module_type := fcase(
  module %in% c("grey", "gold"), "technical",
  default = "biological"
)]
pres_bio <- pres_dt[module_type == "biological"]

me_cols <- setdiff(names(MEs_all), "sample")
r1_meta <- meta[core == r1, .(label, age_kyr = y_bp / 1000)]
r2_meta <- meta[core == r2, .(label, age_kyr = y_bp / 1000)]
r1_me <- MEs_all[sample %in% rownames(expr_by_core[[r1]])]
r2_me <- MEs_all[sample %in% rownames(expr_valid)]

concordance <- rbindlist(lapply(me_cols, function(me) {
  d1 <- merge(r1_meta, r1_me[, .(sample, value = get(me))], by.x = "label", by.y = "sample")
  d2 <- merge(r2_meta, r2_me[, .(sample, value = get(me))], by.x = "label", by.y = "sample")
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

settings_dt <- data.table(
  param = c("soft_power", "deepSplit", "mergeCutHeight", "minModuleSize", "n_perm", "mode"),
  value = c("12", "3", "0.25", "20", as.character(n_perm), run_mode)
)

fwrite(data.table(taxon = names(module_colors), module = as.character(module_colors)),
       assign_file, sep = "\t")
fwrite(MEs_all, file.path(OUT_WGCNA, "module_eigengenes.tsv"), sep = "\t")
fwrite(pres_dt, file.path(OUT_WGCNA, "preservation.tsv"), sep = "\t")
fwrite(pres_bio, file.path(OUT_WGCNA, "preservation_biological.tsv"), sep = "\t")
fwrite(concordance, file.path(OUT_WGCNA, "eigengene_concordance_age_aligned.tsv"), sep = "\t")
fwrite(settings_dt, file.path(OUT_WGCNA, "settings_used.tsv"), sep = "\t")
saveRDS(net, file.path(OUT_WGCNA, "consensus_wgcna.rds"))

message(sprintf(
  "[balancednetwork] balanced WGCNA complete: per_core=%d, modules_non_grey=%d",
  unique(core_n)[1], length(setdiff(unique(module_colors), "grey"))
))
