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
cat("\n")

# 6. Taxon ID Consistency
cat("## Taxon ID Consistency\n")
vst_taxa <- colnames(vst)
mod_taxa <- mods$taxon
meta_taxa <- fread("results/stage1/prokaryotes_taxa_metadata.tsv")$taxon
raw_taxa <- unique(fread(UPSTREAM$tax_damage, select = "subspecies")$subspecies)

cat("- Intersection VST & Modules:", length(intersect(vst_taxa, mod_taxa)), "\n")
cat("- Intersection VST & Metadata:", length(intersect(vst_taxa, meta_taxa)), "\n")
cat("- Intersection VST & Raw (subspecies):", length(intersect(vst_taxa, raw_taxa)), "\n\n")

# 7. Fuzzy Forest Results
if (file.exists("results/importance_fuzzy/ff_model_performance.tsv")) {
  cat("## Fuzzy Forest Results\n")
  ff_perf <- fread("results/importance_fuzzy/ff_model_performance.tsv")
  cat("- Model Accuracy:", ff_perf[Metric == "Accuracy", Value], "\n")
  cat("- Model Kappa:", ff_perf[Metric == "Kappa", Value], "\n\n")
  
  ff_vimp <- fread("results/importance_fuzzy/ff_variable_importance.tsv")
  cat("- Top 5 Fuzzy Forest Drivers:\n\n")
  print(knitr::kable(head(ff_vimp[, .(taxon, variable_importance, module)], 5)))
  cat("\n")
}

# 8. Network Statistics Results
if (file.exists("results/network_stats/network_metrics_summary.tsv")) {
  cat("## Network Statistics Results\n")
  net_stats <- fread("results/network_stats/network_metrics_summary.tsv")
  cat("- Total taxa analyzed:", nrow(net_stats), "\n")
  cat("- Potential keystones identified:", sum(net_stats$is_potential_keystone, na.rm = TRUE), "\n")
  cat("- Importance Type Distribution:\n\n")
  imp_counts <- net_stats[, .N, by = importance_type][order(-N)]
  print(knitr::kable(imp_counts))
  cat("\n")
  
  cat("- Top 5 Keystone Taxa:\n\n")
  keystones <- net_stats[is_potential_keystone == TRUE][order(rank_sum)]
  print(knitr::kable(head(keystones[, .(taxon, module, functional_group, species)], 5)))
  cat("\n")
}

# 9. Integrated Driver Results
if (file.exists("results/importance/integrated_driver_summary.tsv")) {
  cat("## Integrated Driver Results\n")
  drivers <- fread("results/importance/integrated_driver_summary.tsv")
  cat("- Super-Drivers (Statistical + Topological):", sum(drivers$is_super_driver, na.rm = TRUE), "\n")
  cat("- Driver Categories:\n\n")
  cat_counts <- drivers[, .N, by = driver_category][order(-N)]
  print(knitr::kable(cat_counts))
  cat("\n")
  
  cat("- Super-Driver List:\n\n")
  sd_list <- drivers[is_super_driver == TRUE][order(-integrated_score)]
  print(knitr::kable(sd_list[, .(taxon, module, functional_group, genus, variable_importance, pagerank)]))
  cat("\n")
}

sink()
message("Summary written to DATA_SUMMARY.md")
