#!/usr/bin/env Rscript
# 06_emp.R — Encoded Metabolic Potential (EMP) and Sugar/Acid Pathway (SAP)
#
# EMP = Σ_sample [ Σ_taxon( VST_abundance × Σ_module(completeness × |ΔG|) ) ]
#   — thermodynamic-weighted aggregate metabolic capacity per sample
#
# SAP = abundance-weighted mean(sugar_module_completeness / acid_module_completeness)
#   — organic matter decomposition pathway preference
#
# Inputs:  results/stage1/prokaryotes_vst.rds
#          CLASS$prokaryote_function, CLASS$thermo
#          results/stage1/wgcna/module_eigengenes.tsv
#          UPSTREAM$metadata
# Outputs: results/emp/
#   emp_sap_per_sample.tsv      EMP, SAP, EMP_scaled per sample
#   taxon_dg_capacity.tsv       per-taxon thermodynamic capacity
#   taxon_sap.tsv               per-taxon sugar/acid pathway completeness
#   trait_table_with_emp.rds    merged eigengene + EMP table (used by 05_tea_vs_emp.R)
#   emp_sap_summary.tsv         EMP/SAP correlation summary

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── 1. Load data ──────────────────────────────────────────────────────────────

log_msg("Loading data...")

thermo    <- fread(CLASS$thermo)
prok_func <- fread(CLASS$prokaryote_function)
vst       <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
MEs       <- fread(file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv"))
meta      <- fread(UPSTREAM$metadata)

trait_table <- merge(MEs, meta[, .(label, core, y_bp, mis)],
                     by.x = "sample", by.y = "label")
trait_table[, age_kyr := y_bp / 1000]
setnames(trait_table, "mis", "d18O")

log_msg(sprintf("  VST: %d samples × %d taxa", nrow(vst), ncol(vst)))
log_msg(sprintf("  Taxa with function data: %d", nrow(prok_func)))

# ── 2. Parse module completeness per taxon ────────────────────────────────────
# Format in module_completeness_detail: "M00307:1;W00008:0.33;..."

log_msg("Parsing module completeness per taxon...")

parse_modules <- function(detail_str) {
  if (is.na(detail_str) || detail_str == "") return(data.table())
  parts <- strsplit(detail_str, ";", fixed = TRUE)[[1]]
  rbindlist(lapply(parts, function(p) {
    kv <- strsplit(p, ":", fixed = TRUE)[[1]]
    if (length(kv) == 2) data.table(module = kv[1], completeness = as.numeric(kv[2]))
    else data.table()
  }))
}

taxon_modules <- rbindlist(lapply(seq_len(nrow(prok_func)), function(i) {
  mods <- parse_modules(prok_func$module_completeness_detail[i])
  if (nrow(mods) > 0) mods[, taxon := prok_func$taxon[i]]
  mods
}), fill = TRUE)

log_msg(sprintf("  Parsed %d taxon-module pairs", nrow(taxon_modules)))

# ── 3. Map modules to |ΔG| ────────────────────────────────────────────────────

log_msg("Mapping modules to |ΔG| values...")

module_dg <- rbindlist(lapply(seq_len(nrow(thermo)), function(i) {
  mods <- strsplit(thermo$kegg_modules[i], ",", fixed = TRUE)[[1]]
  data.table(module = trimws(mods),
             reaction = thermo$reaction_id[i],
             delta_g  = abs(thermo$delta_G_kJ_mol[i]),
             oxygen_regime = thermo$oxygen_regime[i])
}))

taxon_dg <- merge(taxon_modules, module_dg, by = "module", all.x = TRUE)
# Fallback: taxa whose modules are not in the thermodynamics table get
# emp_default_dg = 100 kJ/mol, a conservative estimate for heterotrophic
# organic matter oxidation (cf. Thauer et al. 1977).
taxon_dg[is.na(delta_g), `:=`(delta_g = PARAMS$emp_default_dg,
                               reaction = "heterotrophy_default")]

# ── 4. Per-taxon thermodynamic capacity ──────────────────────────────────────

taxon_capacity <- taxon_dg[, .(
  total_dg_capacity = sum(completeness * delta_g, na.rm = TRUE),
  n_modules         = .N,
  primary_reaction  = reaction[which.max(completeness * delta_g)]
), by = taxon]

log_msg(sprintf("  Taxon capacity computed for %d taxa", nrow(taxon_capacity)))

# ── 5. EMP per sample ─────────────────────────────────────────────────────────

log_msg("Calculating EMP per sample...")

common_taxa <- intersect(taxon_capacity$taxon, colnames(vst))
log_msg(sprintf("  Taxa in both datasets: %d", length(common_taxa)))

cap_vec <- taxon_capacity[match(common_taxa, taxon), total_dg_capacity]
cap_vec[is.na(cap_vec)] <- PARAMS$emp_default_dg

# CLR values are zero-centred and can be negative. Shift to ≥ 0 before the dot
# product so that low-abundance taxa contribute positively rather than subtracting
# from EMP. The + 0.01 floor avoids exact zeros in the weight vector. This shift
# is applied globally (not per-sample) so relative taxon weights within each
# sample are preserved.
vst_common  <- vst[, common_taxa]
vst_shifted <- vst_common - min(vst_common) + 0.01

emp_dt <- data.table(
  sample    = rownames(vst),
  EMP       = as.vector(vst_shifted %*% cap_vec),
  EMP_scaled = NA_real_
)
emp_dt[, EMP_scaled := EMP / 1e6]

# ── 6. SAP per sample ────────────────────────────────────────────────────────

log_msg("Calculating SAP per sample...")

sugar_modules <- c("M00001", "M00002", "M00003", "M00004", "M00005",
                   "M00006", "M00007", "M00580", "M00579")
acid_modules  <- c("M00307", "M00009", "M00010", "M00011", "M00620", "M00377")

taxon_sugar <- taxon_modules[module %in% sugar_modules,
                              .(sugar_cap = sum(completeness, na.rm = TRUE)), by = taxon]
taxon_acid  <- taxon_modules[module %in% acid_modules,
                              .(acid_cap  = sum(completeness, na.rm = TRUE)), by = taxon]

taxon_sap   <- merge(taxon_sugar, taxon_acid, by = "taxon", all = TRUE)
taxon_sap[is.na(sugar_cap), sugar_cap := 0]
# Floor acid_cap at 0.01 (not 0) to avoid division by zero; completeness scores
# are in [0, 1], so 0.01 is effectively absent without producing Inf.
taxon_sap[is.na(acid_cap),  acid_cap  := 0.01]
taxon_sap[, sap := sugar_cap / acid_cap]

common_sap <- intersect(taxon_sap$taxon, colnames(vst))
sap_vec    <- taxon_sap[match(common_sap, taxon), sap]
sap_vec[is.na(sap_vec)] <- 1

vst_sap         <- vst[, common_sap]
vst_sap_shifted <- vst_sap - min(vst_sap) + 0.01
sap_wt          <- (vst_sap_shifted %*% sap_vec) / rowSums(vst_sap_shifted)
emp_dt[, SAP := as.vector(sap_wt)]

# ── 7. Merge with trait table and correlate with climate ─────────────────────

log_msg("Merging with trait table and computing climate correlations...")

trait_emp <- merge(trait_table, emp_dt, by = "sample", all.x = TRUE)

emp_d18o <- cor.test(trait_emp$EMP, trait_emp$d18O, use = "complete.obs")
sap_d18o <- cor.test(trait_emp$SAP, trait_emp$d18O, use = "complete.obs")
emp_age  <- cor.test(trait_emp$EMP, trait_emp$age_kyr, use = "complete.obs")

log_msg(sprintf("  EMP ~ d18O: r = %.3f, p = %.2e", emp_d18o$estimate, emp_d18o$p.value))
log_msg(sprintf("  SAP ~ d18O: r = %.3f, p = %.2e", sap_d18o$estimate, sap_d18o$p.value))
log_msg(sprintf("  EMP ~ age:  r = %.3f, p = %.2e", emp_age$estimate,  emp_age$p.value))

# ── 8. Save ───────────────────────────────────────────────────────────────────

log_msg("Saving outputs...")

fwrite(emp_dt,         file.path(RESULTS$emp, "emp_sap_per_sample.tsv"), sep = "\t")
fwrite(taxon_capacity, file.path(RESULTS$emp, "taxon_dg_capacity.tsv"),  sep = "\t")
fwrite(taxon_sap,      file.path(RESULTS$emp, "taxon_sap.tsv"),          sep = "\t")
saveRDS(trait_emp,     file.path(RESULTS$emp, "trait_table_with_emp.rds"))

summary_dt <- data.table(
  metric       = c("EMP", "SAP"),
  mean         = c(mean(emp_dt$EMP, na.rm = TRUE),    mean(emp_dt$SAP, na.rm = TRUE)),
  sd           = c(sd(emp_dt$EMP,   na.rm = TRUE),    sd(emp_dt$SAP,   na.rm = TRUE)),
  cor_d18O     = c(emp_d18o$estimate,                 sap_d18o$estimate),
  cor_d18O_p   = c(emp_d18o$p.value,                  sap_d18o$p.value),
  cor_age      = c(emp_age$estimate,                   NA_real_)  # SAP~age not computed
)
fwrite(summary_dt, file.path(RESULTS$emp, "emp_sap_summary.tsv"), sep = "\t")

log_msg("Done. Outputs in ", RESULTS$emp)
