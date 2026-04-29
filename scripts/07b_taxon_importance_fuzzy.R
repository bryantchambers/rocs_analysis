#!/usr/bin/env Rscript
# 07b_taxon_importance_fuzzy.R — Quantify State-Specific Taxon Importance using FuzzyForest
#
# Pipeline:
#   1. Load VST data, WGCNA modules, and HMM states
#   2. Split into Training (70%) and Testing (30%) sets
#   3. Apply Fuzzy Forest (ff) using WGCNA modules as blocks
#   4. Evaluate performance on test set
#   5. Extract Variable Importance (VIMP)
#
# Inputs:  results/stage1/prokaryotes_vst.rds
#          results/stage1/wgcna/module_assignments.tsv
#          results/hmm/hmm_states.tsv
# Outputs: results/importance_fuzzy/
#   ff_variable_importance.tsv  Fuzzy Forest importance ranking
#   ff_model_performance.tsv    Accuracy/Kappa on test set
#   ff_model.rds                fitted fuzzy_forest object

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(WGCNA)
  library(fuzzyforest)
  library(caret)
  library(randomForest)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

IMPORTANCE_OUT <- RESULTS$importance_fuzzy
dir.create(IMPORTANCE_OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────

log_msg("Loading datasets...")

# Using relative paths as verified in inspect_data.R
vst   <- readRDS("results/stage1/prokaryotes_vst.rds")
mods  <- fread("results/stage1/wgcna/module_assignments.tsv")
states <- fread("results/hmm/hmm_states.tsv")

# Align samples
common_samples <- intersect(rownames(vst), states$sample)
log_msg(sprintf("  Aligned on %d samples", length(common_samples)))

vst_sub    <- vst[common_samples, ]
states_sub <- states[match(common_samples, sample), ]

# Ensure taxa in VST and mods match
common_taxa <- intersect(colnames(vst_sub), mods$taxon)
vst_sub <- vst_sub[, common_taxa]
mods_sub <- mods[match(common_taxa, taxon), ]

log_msg(sprintf("  Analyzing %d taxa across %d modules", length(common_taxa), length(unique(mods_sub$module))))

# ── 2. Data Splitting ─────────────────────────────────────────────────────────

log_msg("Splitting data into training and testing sets (70/30)...")

target_state <- factor(make.names(states_sub$label))
train_idx <- createDataPartition(target_state, p = 0.7, list = FALSE)

X_train <- as.data.frame(vst_sub[train_idx, ])
y_train <- target_state[train_idx]

X_test <- as.data.frame(vst_sub[-train_idx, ])
y_test <- target_state[-train_idx]

# ── 3. Fuzzy Forest Configuration ─────────────────────────────────────────────

log_msg("Configuring Fuzzy Forest...")

# Module membership vector
module_membership <- mods_sub$module

# Screening parameters: 
# drop_fraction: proportion of features dropped at each iteration of RFE
# keep_fraction: proportion of features to keep from each module
screen_c <- screen_control(drop_fraction = 0.25, 
                           keep_fraction = 0.05, # Keep top 5% of each module
                           min_ntree = 500,
                           ntree_factor = 1)

# Selection parameters:
# number_selected: total number of features to keep for the final forest
select_c <- select_control(drop_fraction = 0.25,
                           number_selected = 50, # Aim for top 50 drivers
                           min_ntree = 500,
                           ntree_factor = 1)

# ── 4. Fit Fuzzy Forest ───────────────────────────────────────────────────────

log_msg("Fitting Fuzzy Forest (this may take a while)...")

# ff() fits recursive feature elimination within each module block,
# then fits a final random forest on the survivors.
ff_fit <- ff(X = X_train, 
             y = y_train, 
             module_membership = module_membership,
             screen_params = screen_c,
             select_params = select_c,
             final_ntree = 2000,
             num_processors = 1) # Set to 1 for stability in container

log_msg("Fuzzy Forest fit complete.")

# ── 5. Evaluation ─────────────────────────────────────────────────────────────

log_msg("Evaluating on test set...")

# Predict using the final random forest inside the fuzzy_forest object
preds <- predict(ff_fit, new_data = X_test)
conf_mat <- confusionMatrix(preds, y_test)

perf_dt <- data.table(
  Metric = c("Accuracy", "Kappa"),
  Value  = c(conf_mat$overall["Accuracy"], conf_mat$overall["Kappa"])
)
fwrite(perf_dt, file.path(IMPORTANCE_OUT, "ff_model_performance.tsv"), sep = "\t")

log_msg(sprintf("  Test Accuracy: %.3f, Kappa: %.3f", 
                perf_dt$Value[1], perf_dt$Value[2]))

# ── 6. Variable Importance & Results ──────────────────────────────────────────

log_msg("Extracting variable importance...")

# feature_list contains the selected features and their importance
vims <- as.data.table(ff_fit$feature_list)
setnames(vims, "feature_name", "taxon")

# Merge with module assignments for context
vims <- merge(vims, mods_sub, by = "taxon", all.x = TRUE)
setorder(vims, -variable_importance)

fwrite(vims, file.path(IMPORTANCE_OUT, "ff_variable_importance.tsv"), sep = "\t")

# Save model object
saveRDS(ff_fit, file.path(IMPORTANCE_OUT, "ff_model.rds"))

# ── 7. Plotting (Optional Diagnostic) ─────────────────────────────────────────

log_msg("Generating diagnostic plots...")
png(file.path(IMPORTANCE_OUT, "ff_modplot.png"), width = 800, height = 600)
modplot(ff_fit)
dev.off()

log_msg("Done. Results in ", IMPORTANCE_OUT)
