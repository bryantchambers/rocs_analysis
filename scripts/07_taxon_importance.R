#!/usr/bin/env Rscript
# 07_taxon_importance.R — Quantify State-Specific Taxon Importance
#
# Pipeline:
#   1. Intra-modular Hubs: Calculate Module Membership (kME)
#   2. State-Membership Scores: Weight kME by state eigengene fingerprints
#   3. Feature Selection: Random Forest with CV to predict HMM states
#
# Inputs:  results/stage1/prokaryotes_vst.rds
#          results/stage1/wgcna/module_assignments.tsv
#          results/stage1/wgcna/module_eigengenes.tsv
#          results/hmm/hmm_states.tsv
#          results/hmm/state_fingerprints.tsv
# Outputs: results/importance/
#   taxon_kme.tsv               Module membership (kME) for all taxa
#   state_importance_scores.tsv kME-fingerprint weighted importance
#   rf_variable_importance.tsv  Random Forest Gini importance per taxon
#   rf_model_performance.tsv    Accuracy/Kappa on withheld test set
#   rf_model.rds                fitted ranger model object

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
  library(caret)
  library(ranger)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

IMPORTANCE_OUT <- RESULTS$importance
dir.create(IMPORTANCE_OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────

log_msg("Loading datasets...")

vst   <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
mods  <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
MEs   <- fread(file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv"))
states <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))
fingerprints <- fread(file.path(RESULTS$hmm, "state_fingerprints.tsv"))

# Ensure sample alignment across VST, MEs, and states
common_samples <- intersect(rownames(vst), MEs$sample)
common_samples <- intersect(common_samples, states$sample)
log_msg(sprintf("  Aligned on %d samples", length(common_samples)))

vst_sub    <- vst[common_samples, ]
MEs_sub    <- MEs[match(common_samples, sample), ]
states_sub <- states[match(common_samples, sample), ]

# ── 2. Intra-modular Hubs (kME) ────────────────────────────────────────────────

log_msg("Calculating Module Membership (kME)...")

# signedKME correlates taxon abundance with module eigengenes.
# High kME (> 0.8) indicates a "hub" taxon for that module.
me_cols <- grep("^ME", names(MEs_sub), value = TRUE)
datME   <- MEs_sub[, ..me_cols]
kME_mat <- signedKME(vst_sub, datME, outputColumnName = "kME")

kME_dt <- as.data.table(kME_mat, keep.rownames = "taxon")
fwrite(kME_dt, file.path(IMPORTANCE_OUT, "taxon_kme.tsv"), sep = "\t")

# ── 3. State-Membership Scores ────────────────────────────────────────────────

log_msg("Calculating State-Membership Scores...")

# Formula: Score(taxon, state) = Σ_module [ kME(taxon, module) * Fingerprint(state, module) ]
# This identifies taxa that are both central to their module and belong to 
# modules that are highly active in a specific HMM state.

# Prepare fingerprints for multiplication
fp_cols <- grep("^ME", names(fingerprints), value = TRUE)
fp_mat  <- as.matrix(fingerprints[, ..fp_cols])
rownames(fp_mat) <- fingerprints$label

# Match kME columns to fingerprint columns
kme_cols <- paste0("kME", fp_cols)
kME_sub_mat <- as.matrix(kME_dt[, ..kme_cols])
rownames(kME_sub_mat) <- kME_dt$taxon

# Matrix multiplication: [Taxa x Modules] %*% [Modules x States] -> [Taxa x States]
state_scores <- kME_sub_mat %*% t(fp_mat)

state_scores_dt <- as.data.table(state_scores, keep.rownames = "taxon")
fwrite(state_scores_dt, file.path(IMPORTANCE_OUT, "state_importance_scores.tsv"), sep = "\t")

# ── 4. Random Forest Feature Selection ────────────────────────────────────────

log_msg("Fitting Random Forest for state prediction...")

# Target variable: HMM state label
ml_data <- as.data.table(vst_sub)
ml_data[, target_state := factor(states_sub$label)]

# Withholding: Split into training (70%) and testing (30%)
train_idx <- createDataPartition(ml_data$target_state, p = 0.7, list = FALSE)
train_set <- ml_data[train_idx, ]
test_set  <- ml_data[-train_idx, ]

log_msg(sprintf("  Training set: %d samples, Test set: %d samples", 
                nrow(train_set), nrow(test_set)))

# Cross-validation: 5-fold CV on training set
train_control <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  classProbs = TRUE,
  verboseIter = FALSE
)

# Fit Random Forest using ranger (efficient implementation)
# Importance = 'impurity' (Gini index) provides a ranking of taxonomic drivers.
rf_fit <- train(
  target_state ~ ., 
  data = train_set,
  method = "ranger",
  trControl = train_control,
  importance = "impurity",
  num.trees = 500
)

log_msg("  RF Model Training complete.")
print(rf_fit)

# ── 5. Model Validation and Variable Importance ───────────────────────────────

log_msg("Evaluating model on withheld test set...")

preds <- predict(rf_fit, newdata = test_set)
conf_mat <- confusionMatrix(preds, test_set$target_state)

perf_dt <- data.table(
  Metric = c("Accuracy", "Kappa"),
  Value  = c(conf_mat$overall["Accuracy"], conf_mat$overall["Kappa"])
)
fwrite(perf_dt, file.path(IMPORTANCE_OUT, "rf_model_performance.tsv"), sep = "\t")

log_msg(sprintf("  Test Accuracy: %.3f, Kappa: %.3f", 
                perf_dt$Value[1], perf_dt$Value[2]))

# Extract Variable Importance
var_imp <- varImp(rf_fit, scale = FALSE)
imp_dt  <- as.data.table(var_imp$importance, keep.rownames = "taxon")
setorder(imp_dt, -Overall)

# Merge with module assignments for context
imp_dt <- merge(imp_dt, mods, by = "taxon", all.x = TRUE)
fwrite(imp_dt, file.path(IMPORTANCE_OUT, "rf_variable_importance.tsv"), sep = "\t")

# ── 6. Save ───────────────────────────────────────────────────────────────────

log_msg("Saving model and final outputs...")
saveRDS(rf_fit, file.path(IMPORTANCE_OUT, "rf_model.rds"))

log_msg("Done. Results in ", IMPORTANCE_OUT)
