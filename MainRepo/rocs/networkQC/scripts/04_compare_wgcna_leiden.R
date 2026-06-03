#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

OUT_BASE <- here("networkQC", "results")
OUT_TABLE <- file.path(OUT_BASE, "tables")
REPORT <- file.path(OUT_BASE, "NETWORK_QC_REPORT.md")

mods_w <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
mods_l <- fread(file.path(OUT_TABLE, "leiden_module_assignments.tsv"))
sweep <- fread(file.path(OUT_TABLE, "qc_parameter_sweep_summary.tsv"))
basics <- fread(file.path(OUT_TABLE, "qc_basics_summary.tsv"))
pairc <- fread(file.path(OUT_TABLE, "qc_pairwise_core_eigengene_concordance.tsv"))
setnames(pairc, old = names(pairc), new = sub("\\.V1$", "", names(pairc)))
leids <- fread(file.path(OUT_TABLE, "leiden_resolution_summary.tsv"))
leidrun <- fread(file.path(OUT_TABLE, "leiden_run_summary.tsv"))

joined <- merge(mods_w, mods_l, by = "taxon", suffixes = c("_wgcna", "_leiden"))
conf <- joined[, .N, by = .(module_wgcna, module_leiden)][order(-N)]
fwrite(conf, file.path(OUT_TABLE, "wgcna_vs_leiden_confusion.tsv"), sep = "\t")

top_conf <- conf[, head(.SD, 1), by = module_wgcna]
top_conf[, module_size := sum(N), by = module_wgcna]
top_conf[, purity := N / module_size]
fwrite(top_conf, file.path(OUT_TABLE, "wgcna_vs_leiden_module_purity.tsv"), sep = "\t")

sink(REPORT)
cat("# Network QC Report\n\n")
cat("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("- Scope: WGCNA QC baseline, parameter sweep, Leiden comparison\n\n")

cat("## 1. Baseline WGCNA QC\n\n")
cat("Key metrics:\n\n")
cat("|metric|value|\n|---|---:|\n")
for (i in seq_len(nrow(basics))) cat(sprintf("|%s|%s|\n", basics$metric[i], as.character(basics$value[i])))
cat("\n")

cat("## 2. Parameter Sweep Summary\n\n")
ok <- sweep[status == "ok"]
cat(sprintf("- Total combinations tested: %d\n", nrow(sweep)))
cat(sprintf("- Successful fits: %d\n", nrow(ok)))
if (nrow(ok) > 0) {
  best <- ok[order(grey_pct, -non_grey_modules)][1]
  cat(sprintf("- Lowest grey_pct candidate: power=%d, deepSplit=%d, mergeCutHeight=%.2f, minModuleSize=%d, grey_pct=%.2f\n\n",
              best$power, best$deepSplit, best$mergeCutHeight, best$minModuleSize, best$grey_pct))
}

cat("## 3. Cross-core Eigengene Concordance\n\n")
pair_sum <- pairc[is.finite(pearson_r), .(
  pearson_mean = mean(pearson_r),
  spearman_mean = mean(spearman_rho),
  rmse_mean = mean(rmse)
), by = .(core_a, core_b)][order(-pearson_mean)]
cat("|core_a|core_b|pearson_mean|spearman_mean|rmse_mean|\n|---|---|---:|---:|---:|\n")
for (i in seq_len(nrow(pair_sum))) {
  cat(sprintf("|%s|%s|%.3f|%.3f|%.3f|\n",
              pair_sum$core_a[i], pair_sum$core_b[i], pair_sum$pearson_mean[i], pair_sum$spearman_mean[i], pair_sum$rmse_mean[i]))
}
cat("\n")

cat("## 4. Leiden Comparison\n\n")
cat("Resolution scan:\n\n")
cat("|resolution|n_modules|module_size_median|largest_module|\n|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(leids))) {
  cat(sprintf("|%.2f|%d|%.1f|%d|\n", leids$resolution[i], leids$n_modules[i], leids$module_size_median[i], leids$largest_module[i]))
}
cat("\nSelected Leiden run:\n\n")
cat("|metric|value|\n|---|---|\n")
for (i in seq_len(nrow(leidrun))) cat(sprintf("|%s|%s|\n", leidrun$metric[i], as.character(leidrun$value[i])))
cat("\n")

cat("WGCNA-to-Leiden dominant mapping purity:\n\n")
cat("|module_wgcna|module_leiden|N|purity|\n|---|---|---:|---:|\n")
for (i in seq_len(min(20, nrow(top_conf[order(-purity)])))) {
  r <- top_conf[order(-purity)][i]
  cat(sprintf("|%s|%s|%d|%.3f|\n", r$module_wgcna, r$module_leiden, r$N, r$purity))
}
cat("\n")

cat("## 5. Next Actions\n\n")
cat("1. Use sweep results to choose 2-3 candidate WGCNA settings with lower grey burden.\n")
cat("2. Re-run stability (`02b`) for each candidate and compare Jaccard + concordance.\n")
cat("3. Compare downstream outcomes (07b/09/10/11/12) under WGCNA-best vs Leiden-best module sets.\n")
sink()

message("[networkQC] report written: ", REPORT)
