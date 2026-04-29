#!/usr/bin/env Rscript
# scripts/inspect_data.R - Summarize key data for taxon importance analysis

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

output_file <- here("DATA_SUMMARY.md")

sink(output_file)
cat("# Data Summary for Taxon Importance Analysis\n\n")

# 1. VST Data
cat("## Prokaryote VST Data\n")
vst <- readRDS("results/stage1/prokaryotes_vst.rds")
cat("- Dimensions:", paste(dim(vst), collapse = " x "), "(Samples x Taxa)\n")
cat("- Range of values:", paste(range(vst, na.rm = TRUE), collapse = " to "), "\n")
cat("- Any NAs:", any(is.na(vst)), "\n\n")

# 2. Module Assignments
cat("## WGCNA Module Assignments\n")
mods <- fread("results/stage1/wgcna/module_assignments.tsv")
cat("- Total taxa assigned:", nrow(mods), "\n")
cat("- Module Distribution:\n\n")
mod_counts <- mods[, .N, by = module][order(-N)]
print(knitr::kable(mod_counts))
cat("\n")

# 3. HMM States
cat("## HMM States\n")
states <- fread("results/hmm/hmm_states.tsv")
cat("- Total samples with state labels:", nrow(states), "\n")
cat("- State Label Distribution:\n\n")
state_counts <- states[, .N, by = label][order(label)]
print(knitr::kable(state_counts))
cat("\n")

# 4. Sample Alignment
cat("## Sample Alignment\n")
common_samples <- intersect(rownames(vst), states$sample)
cat("- Samples in both VST and HMM:", length(common_samples), "\n")
cat("- Samples missing from VST:", length(setdiff(states$sample, rownames(vst))), "\n")
cat("- Samples missing from HMM:", length(setdiff(rownames(vst), states$sample)), "\n\n")

# 5. FuzzyForest specific check
cat("## FuzzyForest Requirements\n")
cat("- Features with zero variance:\n")
zv <- which(apply(vst, 2, var) == 0)
if(length(zv) > 0) {
  cat("  - Count:", length(zv), "\n")
} else {
  cat("  - None\n")
}

sink()
message("Summary written to DATA_SUMMARY.md")
