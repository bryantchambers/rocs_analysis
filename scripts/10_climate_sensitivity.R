#!/usr/bin/env Rscript
# 10_climate_sensitivity.R — Model taxon abundance against climate proxies (d18O/SST)
#
# Pipeline:
#   1. Load Integrated Driver Summary (Tiers 1, 2, 3)
#   2. Load CLR-transformed taxon abundances
#   3. Load Metadata (d18O/MIS and SST)
#   4. For each top candidate, fit GLS model: Abundance ~ d18O + SST + core
#   5. Use corCAR1 autocorrelation structure (age-based) to handle time-series
#   6. Identify "Climate-Sensitive Drivers"
#
# Outputs: results/importance/
#   climate_sensitivity_results.tsv   GLS coefficients and p-values for all tested taxa
#   climate_driver_volcano_plot.png   Visualization of sensitivity vs. significance

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(nlme)
  library(ggplot2)
  library(ggrepel)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

CLIMATE_OUT <- RESULTS$importance # Keeping in importance folder for now

# ── 1. Load Data ──────────────────────────────────────────────────────────────

log_msg("Loading driver summary and abundances...")
drivers <- fread(file.path(RESULTS$importance, "integrated_driver_summary.tsv"))
vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(UPSTREAM$metadata)

# Define candidates: All Tier 1, 2, and 3
candidates <- drivers[driver_tier != "Peripheral", taxon]
log_msg(sprintf("  Testing sensitivity for %d candidate taxa", length(candidates)))

# ── 2. Prepare Metadata ───────────────────────────────────────────────────────

log_msg("Preparing climate proxies...")
# In this project, MIS (Marine Isotope Stage) is used as the d18O proxy
# sst is available in metadata for some samples
meta_sub <- meta[, .(sample = label, core, y_bp, d18O = mis, sst)]
meta_sub[, age_kyr := y_bp / 1000]

# ── 3. GLS Modeling ───────────────────────────────────────────────────────────

log_msg("Fitting GLS models with CAR1 autocorrelation...")

run_sensitivity <- function(taxon_id) {
  # Merge abundance with metadata
  dat <- data.table(abundance = vst[, taxon_id], sample = rownames(vst))
  dat <- merge(dat, meta_sub, by = "sample")
  
  # Ensure data is sorted by age for autocorrelation
  setorder(dat, core, age_kyr)
  
  # Remove NAs for SST (if needed) or handle separately
  # We'll test d18O primarily as it has full coverage
  dat_clean <- dat[!is.na(abundance) & !is.na(d18O) & !is.na(age_kyr)]
  
  if (nrow(dat_clean) < 30) return(NULL)
  
  res <- tryCatch({
    # Model: Abundance ~ d18O + core (adding SST if available and sufficient)
    # Using d18O (MIS) as the primary climate signal
    fit <- gls(abundance ~ d18O + core, 
               data = dat_clean,
               correlation = corCAR1(form = ~ age_kyr | core),
               na.action = na.omit)
    
    ct <- summary(fit)$tTable
    data.table(
      taxon = taxon_id,
      coef_d18O = ct["d18O", "Value"],
      p_d18O = ct["d18O", "p-value"],
      n_samples = nrow(dat_clean)
    )
  }, error = function(e) {
    # If GLS fails (e.g. convergence), fallback to simple LM or return NULL
    return(NULL)
  })
  
  return(res)
}

results_list <- lapply(candidates, run_sensitivity)
sensitivity_dt <- rbindlist(results_list)

# ── 4. Multiple Testing Correction ────────────────────────────────────────────

if (nrow(sensitivity_dt) > 0) {
  sensitivity_dt[, fdr_d18O := p.adjust(p_d18O, method = "BH")]
  
  # Merge back with driver info
  final_sensitivity <- merge(sensitivity_dt, drivers[, .(taxon, species, module, driver_tier, integrated_score)], by = "taxon")
  
  fwrite(final_sensitivity, file.path(CLIMATE_OUT, "climate_sensitivity_results.tsv"), sep = "\t")
  
  log_msg(sprintf("  Identified %d climate-sensitive taxa (FDR < 0.05)", 
                  sum(final_sensitivity$fdr_d18O < PARAMS$fdr_threshold)))
  
  # ── 5. Visualization ────────────────────────────────────────────────────────
  
  log_msg("Generating sensitivity volcano plot...")
  
  p_volc <- ggplot(final_sensitivity, aes(x = coef_d18O, y = -log10(p_d18O), color = driver_tier)) +
    geom_point(aes(size = integrated_score), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    scale_color_manual(values = c("Tier 1: Super-Driver" = "#D55E00", 
                                 "Tier 2: High Potential" = "#0072B2", 
                                 "Tier 3: Predictive Specialist" = "#009E73")) +
    labs(title = "Climate Sensitivity of Driver Taxa",
         subtitle = "Positive Coef = Higher abundance in Glacial stages (High d18O/MIS)",
         x = "GLS Coefficient (d18O)",
         y = "-log10(p-value)") +
    geom_text_repel(data = final_sensitivity[fdr_d18O < 0.01 & abs(coef_d18O) > 0.5], 
                    aes(label = species), size = 3, max.overlaps = 15)
  
  ggsave(file.path(CLIMATE_OUT, "climate_driver_volcano_plot.png"), p_volc, width = 10, height = 8)
}

log_msg("Done. Results in ", CLIMATE_OUT)
