#!/usr/bin/env Rscript
# 07_tea_vs_emp.R — TEA indices for ROCS data vs EMP/SAP comparison
#
# Implements the same TEA (Terminal Electron Acceptor) indices used in the
# Mediterranean sapropel study (14-tea-indices.R) and compares them to the
# existing EMP/SAP metrics computed in 06_emp.R.
#
# ── Index definitions (following 14-tea-indices.R) ───────────────────────────
#
#   OAP   Oxidoreductase Acceptor Potential: abundance-weighted mean E0'
#         Uses module completeness × disc × E0' (vs KO-based approach in sapropel)
#
#   MGI   Methanogenesis Index: log(methanogenesis module completeness)
#   MFI   Methane Filter Index: log(methane oxidation / methanogenesis completeness)
#   DCI   Denitrification Completeness Index: log(M00530 / M00529 completeness)
#         M00530 carries nosZ (N2O reductase, final step); M00529 carries narG
#   MII   Microoxic Interface Index: log(bd-type / aa3-type oxidase completeness)
#         bd-type = microoxic (cydA); aa3-type = aerobic (coxAC)
#   SRPI  Sulfate Reduction Potential Index: M00596 completeness (sulfate red.)
#
# ── Deviation from 14-tea-indices.R ─────────────────────────────────────────
#
#   The sapropel script derives indices from raw KO counts (enzyme_hits_in_module).
#   ROCS data has per-taxon module completeness but not individual KO counts in
#   the classification file. We therefore approximate using stepwise module
#   completeness from UPSTREAM$kegg_mods (kegg-modules-summary*.tsv.gz), which
#   provides the same `enzyme_hits_in_module` field for any taxon whose genome
#   is in the functional database. Taxa not matched fall back to the classification
#   file completeness scores.
#
# ── Comparison ───────────────────────────────────────────────────────────────
#
#   After computing both sets of indices we assess:
#   1. Correlation of each TEA index with EMP and SAP
#   2. Correlation with δ¹⁸O (climate signal)
#   3. Whether TEA ratio indices capture variance beyond the scalar EMP
#
# Inputs:  UPSTREAM$kegg_mods, UPSTREAM$tax_damage, UPSTREAM$metadata
#          CLASS$prokaryote_function
#          results/stage1/prokaryotes_vst.rds
#          results/emp/emp_sap_per_sample.tsv
# Outputs: results/tea/
#   tea_indices_per_sample.tsv   OAP, DCI, MII, SRPI, MGI, MFI per sample
#   oap_per_taxon.tsv            per-taxon OAP score + oxy_class assignment
#   tea_vs_emp_correlations.tsv  pairwise Pearson correlations among all indices
#   tea_climate_models.tsv       GLS results (TEA index ~ d18O, corCAR1 autocorr)

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(nlme)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

OAP_MODS <- rbindlist(lapply(names(TEA$oap_modules), function(m) {
  x <- TEA$oap_modules[[m]]
  data.table(module = m, class = x$class, eo_mv = x$eo_mv, disc = x$disc)
}))

# ── 1. Load and parse KEGG module completeness ────────────────────────────────
# Primary source: UPSTREAM$kegg_mods (has enzyme_hits_in_module for KO-level)

log_msg("Loading KEGG module data from upstream pipeline...")

all_mods_needed <- unique(c(OAP_MODS$module,
                              "M00174",  # methane monooxygenase (MFI numerator)
                              "M00596",  # sulfate reduction
                              "M00567", "M00357", "M00356", "M00563",  # methanogenesis
                              "M00529", "M00530",  # denitrification
                              "M00155", "M00156", "M00153", "M00154",  # aerobic ox.
                              "M00416", "M00417", "M00973"))

kegg_raw <- fread(
  cmd = sprintf("zcat '%s'", UPSTREAM$kegg_mods),
  select = c("module", "genome_name", "stepwise_module_completeness", "enzyme_hits_in_module"),
  showProgress = FALSE
)
setnames(kegg_raw, c("module", "reference", "completeness", "enzyme_hits"))

kegg_filt <- kegg_raw[module %in% all_mods_needed]
log_msg(sprintf("  Loaded %d genome-module rows for relevant modules", nrow(kegg_filt)))

# ── 2. Map genome names to ROCS taxon IDs ─────────────────────────────────────
# ROCS taxa: S__GCA_001940725.1 → genome_id = GCA_001940725.1
# KEGG genome_name: numeric MAG IDs (3300007352_16) or sometimes GCA/GCF strings
# The upstream pipeline uses numeric IDs for MAGs and GCF for NCBI isolates.
# The classification file genome_id column holds GCA for NCBI isolates.

log_msg("Mapping genome names to ROCS taxon IDs...")

prok_func <- fread(CLASS$prokaryote_function,
                   select = c("taxon", "genome_id", "module_completeness_detail"))

# Build a genome_name → taxon lookup
# For NCBI isolates: GCA_*.1 → GCF_*.1 (same numeric body, different prefix)
# Many entries are already direct matches (MAG numeric IDs)

prok_func[, ref_clean := genome_id]   # direct match for MAG IDs

# GCA → GCF conversion (NCBI RefSeq uses GCF)
prok_func[grepl("^GCA_", genome_id),
          gcf_id := sub("^GCA_", "GCF_", genome_id)]

# Build lookup: try direct match first, then GCF conversion
ref_lookup <- rbind(
  prok_func[, .(reference = ref_clean, taxon)],
  prok_func[!is.na(gcf_id), .(reference = gcf_id, taxon)],
  fill = TRUE
)
ref_lookup <- unique(ref_lookup[!is.na(reference)])

kegg_taxa <- merge(kegg_filt, ref_lookup, by = "reference", all.x = TRUE)
matched   <- kegg_taxa[!is.na(taxon)]
log_msg(sprintf("  Matched %d / %d genome-module rows to ROCS taxa",
                nrow(matched), nrow(kegg_filt)))

# Fallback: for unmatched taxa, parse module completeness from classification file
log_msg("Adding fallback module completeness from classification file...")

parse_modules_dt <- function(func_dt) {
  rbindlist(lapply(seq_len(nrow(func_dt)), function(i) {
    s <- func_dt$module_completeness_detail[i]
    if (is.na(s) || s == "") return(data.table())
    parts <- strsplit(s, ";", fixed = TRUE)[[1]]
    mods <- rbindlist(lapply(parts, function(p) {
      kv <- strsplit(p, ":", fixed = TRUE)[[1]]
      if (length(kv) == 2) data.table(module = kv[1], completeness = as.numeric(kv[2]))
      else data.table()
    }))
    if (nrow(mods) > 0) mods[, taxon := func_dt$taxon[i]]
    mods
  }), fill = TRUE)
}

fallback_taxa <- setdiff(prok_func$taxon, matched$taxon)
if (length(fallback_taxa) > 0) {
  fb_dt   <- parse_modules_dt(prok_func[taxon %in% fallback_taxa])
  fb_filt <- fb_dt[module %in% all_mods_needed]
  fb_filt[, `:=`(reference = taxon, enzyme_hits = NA_character_)]
  setcolorder(fb_filt, c("module", "reference", "completeness", "enzyme_hits", "taxon"))
  all_taxon_mods <- rbind(matched[, .(module, taxon, completeness, enzyme_hits)],
                           fb_filt[, .(module, taxon, completeness, enzyme_hits)])
  log_msg(sprintf("  Added %d fallback rows for %d unmatched taxa",
                  nrow(fb_filt), length(fallback_taxa)))
} else {
  all_taxon_mods <- matched[, .(module, taxon, completeness, enzyme_hits)]
}

log_msg(sprintf("  Total: %d taxon-module rows for %d unique taxa",
                nrow(all_taxon_mods), uniqueN(all_taxon_mods$taxon)))

# ── 3. Extract KO counts from enzyme_hits ─────────────────────────────────────
# enzyme_hits_in_module is comma-separated KO list (e.g. "K00399,K00401,K00402")

log_msg("Extracting KO presence from enzyme_hits_in_module...")

ko_from_hits <- all_taxon_mods[!is.na(enzyme_hits) & enzyme_hits != "",
                                .(taxon, module, enzyme_hits)][
  , .(ko = unlist(strsplit(enzyme_hits, ",", fixed = TRUE))), by = .(taxon, module)
][ko %in% TEA$target_kos]

# Count KOs per taxon (presence/absence within genome)
ko_wide <- dcast(ko_from_hits[, .N, by = .(taxon, ko)],
                 taxon ~ ko, value.var = "N", fill = 0)

log_msg(sprintf("  Taxa with target KOs: %d", nrow(ko_wide)))

# ── 4. Load abundance data ────────────────────────────────────────────────────

log_msg("Loading abundance and sample metadata...")

vst  <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(UPSTREAM$metadata)
meta <- meta[label %in% rownames(vst)]

# VST → relative abundance proxy (shift to ≥ 0, then row-normalise)
vst_pos <- vst - min(vst)
rel_ab   <- vst_pos / rowSums(vst_pos)   # row-normalised (sums to 1 per sample)

taxa_in_vst <- colnames(vst)

# ── 5. OAP per taxon ──────────────────────────────────────────────────────────

log_msg("Computing OAP per taxon...")

# Use module completeness for TEA class membership
oap_taxon <- merge(
  all_taxon_mods[module %in% OAP_MODS$module],
  OAP_MODS, by = "module"
# oap_completeness_threshold (= 0.30): exclude taxa where the TEA pathway is
# < 30% complete — partial pathways should not contribute to OAP classification.
)[completeness >= PARAMS$oap_completeness_threshold,
  .(eff = max(completeness * disc)), by = .(taxon, class, eo_mv)]

ref_oap <- oap_taxon[, .(
  oap_v3      = sum(eff * eo_mv) / sum(eff),
  n_classes   = uniqueN(class),
  dominant_class = class[which.max(eff * eo_mv)],
  o2_eff   = sum(eff[class == "O2"],  na.rm = TRUE),
  so4_eff  = sum(eff[class == "SO4"], na.rm = TRUE),
  ch4_eff  = sum(eff[class == "CH4"], na.rm = TRUE),
  no3_eff  = sum(eff[class == "NO3"], na.rm = TRUE)
), by = taxon]

ref_oap[, oxy_class := fcase(
  o2_eff  > 0 & so4_eff == 0 & ch4_eff == 0, "aerobic",
  so4_eff > 0 & o2_eff  == 0,                "sulfate_reducer",
  ch4_eff > 0 & o2_eff  == 0,                "methanogen",
  no3_eff > 0 & o2_eff  == 0,                "denitrifier",
  o2_eff  > 0 & (so4_eff > 0 | ch4_eff > 0), "facultative",
  default = "no_call"
)]

log_msg(sprintf("  Taxa with OAP: %d", nrow(ref_oap)))
log_msg(sprintf("  oxy_class: %s",
                paste(names(table(ref_oap$oxy_class)),
                      as.integer(table(ref_oap$oxy_class)), sep = "=", collapse = ", ")))

# ── 6. Module completeness summaries per taxon for ratio indices ──────────────

log_msg("Summarising module completeness for ratio indices...")

# Best (max) completeness per taxon per module
best_comp <- all_taxon_mods[, .(best_comp = max(completeness, na.rm = TRUE)), by = .(taxon, module)]
comp_wide  <- dcast(best_comp, taxon ~ module, value.var = "best_comp", fill = 0)

# Helper: retrieve column or zeros
get_col <- function(dt, col) if (col %in% names(dt)) dt[[col]] else rep(0, nrow(dt))

# ── 7. Compute TEA indices per sample ─────────────────────────────────────────

log_msg("Computing TEA indices per sample...")

eps_safe <- function(x) pmax(x, min(x[x > 0], na.rm = TRUE) / 2, 1e-6, na.rm = TRUE)

results <- lapply(rownames(vst), function(samp) {
  ra <- rel_ab[samp, ]   # named relative abundance vector for this sample

  # Align to taxa in comp_wide
  taxa_match <- intersect(comp_wide$taxon, names(ra))
  if (length(taxa_match) == 0) return(NULL)

  ra_m     <- ra[taxa_match]
  comp_sub <- comp_wide[taxon %in% taxa_match]
  setorder(comp_sub, taxon)
  ra_m     <- ra_m[comp_sub$taxon]

  # OAP: abundance-weighted mean E0'
  oap_join <- ref_oap[taxon %in% taxa_match]
  if (nrow(oap_join) > 0) {
    oap_join[, ra_i := ra_m[taxon]]
    OAP <- weighted.mean(oap_join$oap_v3, oap_join$ra_i, na.rm = TRUE)
  } else {
    OAP <- NA_real_
  }

  # KO-based ratio indices (preferred path — uses enzyme_hits_in_module counts).
  # e() returns half the smallest positive value observed, providing a
  # data-adaptive pseudocount that avoids log(0) without fixed-floor bias.
  ko_join <- ko_wide[taxon %in% taxa_match]
  if (nrow(ko_join) > 0) {
    ko_join[, ra_i := ra_m[taxon]]
    wt_sum <- function(ko) {
      if (!ko %in% names(ko_join)) return(0)
      sum(ko_join[[ko]] * ko_join$ra_i, na.rm = TRUE)
    }
    Mcr  <- wt_sum("K00399") + wt_sum("K00401") + wt_sum("K00402")
    Pmo  <- wt_sum("K14080")
    Nar  <- wt_sum("K00370")
    Nor  <- wt_sum("K02305")
    Nos  <- wt_sum("K00376")
    Cyd  <- wt_sum("K00425")
    Cox  <- wt_sum("K02274") + wt_sum("K02276")
    Sred <- wt_sum("K00394") + wt_sum("K00395") + wt_sum("K00958") +
            wt_sum("K11180") + wt_sum("K11181")
    M00596_wt <- sum(get_col(comp_sub, "M00596") * ra_m, na.rm = TRUE)

    e <- function(...) {
      vals <- c(...); vals <- vals[vals > 0]
      if (length(vals) == 0) 1e-6 else min(vals) / 2
    }

    MGI  <- log(Mcr  + e(Mcr))
    MFI  <- log((Pmo + e(Pmo)) / (Mcr + e(Mcr)))
    DCI  <- log((Nos + e(Nos)) / (Nar + Nor + e(Nar)))
    MII  <- log((Cyd + e(Cyd)) / (Cox + e(Cox)))
    # SRPI combines module completeness (M00596) and KO evidence (Sred) as the
    # arithmetic mean of two log terms — averaging pathway-level and gene-level
    # evidence for sulfate reduction.
    SRPI <- (log(M00596_wt + 0.01) + log(Sred + e(Sred))) / 2
  } else {
    # Fallback path: no enzyme_hits available; approximate from module completeness.
    # Note: epsilon strategy differs from KO path — fixed 1e-6 floor instead of
    # data-adaptive half-minimum, because completeness values are already bounded
    # in [0, 1] and a fixed floor is numerically stable. SRPI uses a single log
    # term (module completeness only) since KO evidence is unavailable.
    wm <- function(mod) sum(get_col(comp_sub, mod) * ra_m, na.rm = TRUE)
    meth <- wm("M00567") + wm("M00357") + wm("M00356") + wm("M00563")
    ox   <- wm("M00174")  # methane monooxygenase
    no3a <- wm("M00529")  # narG step
    no3b <- wm("M00530")  # nosZ step
    aer  <- wm("M00155") + wm("M00154") + wm("M00416") + wm("M00417")
    mox  <- wm("M00156") + wm("M00153")  # microoxic bd-type
    so4  <- wm("M00596")

    e <- function(x) max(x, 1e-6)
    MGI  <- log(meth + e(meth))
    MFI  <- log((ox + 1e-6) / (meth + 1e-6))
    DCI  <- log((no3b + 1e-6) / (no3a + 1e-6))
    MII  <- log((mox  + 1e-6) / (aer  + 1e-6))
    SRPI <- log(so4 + 0.01)
  }

  data.table(sample = samp, OAP = OAP, MGI = MGI, MFI = MFI,
             DCI = DCI, MII = MII, SRPI = SRPI)
})

tea_dt <- rbindlist(results[!sapply(results, is.null)])
log_msg(sprintf("  TEA indices computed for %d samples", nrow(tea_dt)))

# ── 8. Merge with metadata and EMP ───────────────────────────────────────────

log_msg("Merging with metadata and EMP...")

emp_dt <- fread(file.path(RESULTS$emp, "emp_sap_per_sample.tsv"))

full_dt <- Reduce(function(a, b) merge(a, b, by = "sample", all.x = TRUE),
                  list(tea_dt,
                       emp_dt[, .(sample, EMP, EMP_scaled, SAP)],
                       meta[, .(label, core, y_bp, mis, sst)][ #edited from antionios original see MATEU FOR EDITS!!!
                         , sample := label][, label := NULL]))
full_dt[, age_kyr := y_bp / 1000]
setnames(full_dt, "mis", "d18O")

# ── 9. Climate correlations (GLS with CAR1 autocorrelation) ──────────────────

log_msg("Computing climate correlations (GLS)...")

index_cols <- c("OAP", "MGI", "MFI", "DCI", "MII", "SRPI", "EMP", "SAP")
index_cols <- intersect(index_cols, names(full_dt))

gls_results <- rbindlist(lapply(index_cols, function(idx) {
  sub_dt <- full_dt[is.finite(get(idx)) & !is.na(d18O) & !is.na(age_kyr)]
  if (nrow(sub_dt) < 20) return(NULL)
  tryCatch({
    fit <- gls(as.formula(sprintf("%s ~ d18O + core", idx)),
               data       = sub_dt,
               correlation = corCAR1(form = ~ age_kyr | core),
               na.action   = na.omit)
    ct <- summary(fit)$tTable
    data.table(index = idx,
               coef_d18O = ct["d18O", "Value"],
               se_d18O   = ct["d18O", "Std.Error"],
               t_d18O    = ct["d18O", "t-value"],
               p_d18O    = ct["d18O", "p-value"],
               n         = nrow(sub_dt))
  }, error = function(e) {
    # GLS failed (e.g. CAR1 optimisation did not converge). Fall back to plain
    # Pearson r as a descriptive measure only — p-value is NOT corrected for
    # autocorrelation and should not be used for inference. Check tea_climate_models.tsv
    # for NA p-values to identify which indices triggered this fallback.
    warning(sprintf("GLS failed for %s: %s — using Pearson fallback", idx, conditionMessage(e)))
    r <- cor(sub_dt[[idx]], sub_dt$d18O, use = "complete.obs")
    data.table(index = idx, coef_d18O = r, se_d18O = NA_real_, t_d18O = NA_real_,
               p_d18O = NA_real_, n = nrow(sub_dt))
  })
}))

if (!is.null(gls_results) && nrow(gls_results) > 0) {
  gls_results[, fdr := p.adjust(p_d18O, method = "BH")]
  log_msg("  Climate associations:")
  for (i in seq_len(nrow(gls_results))) {
    log_msg(sprintf("    %-6s β=%.3f  p=%.3e  FDR=%.3f",
                    gls_results$index[i], gls_results$coef_d18O[i],
                    gls_results$p_d18O[i], gls_results$fdr[i]))
  }
}

# ── 10. Pairwise correlations among all indices ────────────────────────────────

log_msg("Computing pairwise correlations among indices...")

corr_mat <- full_dt[, ..index_cols][, lapply(.SD, as.numeric)]
r_mat    <- cor(corr_mat, use = "pairwise.complete.obs")

corr_dt <- as.data.table(as.table(r_mat), keep.rownames = FALSE)
setnames(corr_dt, c("index1", "index2", "r"))
corr_dt <- corr_dt[index1 < index2]  # upper triangle only

# ── 11. Save ──────────────────────────────────────────────────────────────────

log_msg("Saving outputs...")

fwrite(tea_dt,    file.path(RESULTS$tea, "tea_indices_per_sample.tsv"), sep = "\t")
fwrite(ref_oap,   file.path(RESULTS$tea, "oap_per_taxon.tsv"),          sep = "\t")
fwrite(corr_dt,   file.path(RESULTS$tea, "tea_vs_emp_correlations.tsv"), sep = "\t")
if (!is.null(gls_results))
  fwrite(gls_results, file.path(RESULTS$tea, "tea_climate_models.tsv"), sep = "\t")

log_msg("")
log_msg("=== SUMMARY ===")
log_msg(sprintf("  Samples with TEA indices: %d", nrow(tea_dt)))
if (!is.null(gls_results)) {
  log_msg("  Indices significantly associated with d18O (FDR < 0.05):")
  sig <- gls_results[!is.na(fdr) & fdr < PARAMS$fdr_threshold, index]
  log_msg(sprintf("    %s", if (length(sig) > 0) paste(sig, collapse = ", ") else "none"))
}
log_msg("")
log_msg("  Top correlations with EMP:")
emp_corrs <- corr_dt[index1 == "EMP" | index2 == "EMP"][order(-abs(r))]
for (i in seq_len(min(5, nrow(emp_corrs))))
  log_msg(sprintf("    EMP ~ %-6s r = %.3f",
                  ifelse(emp_corrs$index1[i] == "EMP", emp_corrs$index2[i], emp_corrs$index1[i]),
                  emp_corrs$r[i]))

log_msg("Done. Outputs in ", RESULTS$tea)
