#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)

ensure_dirs(c(DIRS$main_tea))

INDEX_COLS <- c("OAP", "MGI", "MFI", "DCI", "MII", "SRPI")
BOUNDARY_VALS <- c(-13.8155105579643, 0)
EPS <- 1e-8
MIN_N_MODULE_TEST <- 20L
MIN_N_STATE_MODULE_TEST <- 8L
MIN_N_STATE_EFFECT_PER_STATE <- 5L

required_inputs <- c(
  INPUTS$tax_damage,
  INPUTS$kegg_mods,
  INPUTS$prokaryote_function,
  file.path(DIRS$main, "sample_metadata_main.tsv"),
  file.path(DIRS$main, "taxa_after_filter.tsv"),
  file.path(DIRS$main, "module_assignments.tsv"),
  file.path(DIRS$main, "hmm_states_main.tsv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required TEA inputs: ", paste(missing_inputs, collapse = "; "))
}

get_col <- function(dt, col) if (col %in% names(dt)) dt[[col]] else rep(0, nrow(dt))

half_min_positive <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) 1e-6 else min(x) / 2
}

parse_modules_dt <- function(func_dt) {
  out <- rbindlist(lapply(seq_len(nrow(func_dt)), function(i) {
    s <- func_dt$module_completeness_detail[i]
    if (is.na(s) || s == "") return(data.table())
    parts <- strsplit(s, ";", fixed = TRUE)[[1]]
    mods <- rbindlist(lapply(parts, function(p) {
      kv <- strsplit(p, ":", fixed = TRUE)[[1]]
      if (length(kv) != 2) return(data.table())
      data.table(module = kv[1], completeness = suppressWarnings(as.numeric(kv[2])))
    }), fill = TRUE)
    if (nrow(mods) > 0) mods[, taxon := func_dt$taxon[i]]
    mods
  }), fill = TRUE)
  out[is.finite(completeness)]
}

compute_indices_ko <- function(ko_join, comp_sub, w_comp) {
  wt_sum <- function(ko) {
    if (!ko %in% names(ko_join)) return(0)
    sum(ko_join[[ko]] * ko_join$w, na.rm = TRUE)
  }

  Mcr <- wt_sum("K00399") + wt_sum("K00401") + wt_sum("K00402")
  Pmo <- wt_sum("K14080")
  Nar <- wt_sum("K00370")
  Nor <- wt_sum("K02305")
  Nos <- wt_sum("K00376")
  Cyd <- wt_sum("K00425")
  Cox <- wt_sum("K02274") + wt_sum("K02276")
  Sred <- wt_sum("K00394") + wt_sum("K00395") + wt_sum("K00958") + wt_sum("K11180") + wt_sum("K11181")
  M00596_wt <- sum(get_col(comp_sub, "M00596") * w_comp, na.rm = TRUE)

  e_Mcr <- half_min_positive(Mcr)
  e_Pmo <- half_min_positive(Pmo)
  e_Nar <- half_min_positive(Nar)
  e_Nos <- half_min_positive(Nos)
  e_Cyd <- half_min_positive(Cyd)
  e_Cox <- half_min_positive(Cox)
  e_Sred <- half_min_positive(Sred)

  data.table(
    MGI = log(Mcr + e_Mcr),
    MFI = log((Pmo + e_Pmo) / (Mcr + e_Mcr)),
    DCI = log((Nos + e_Nos) / (Nar + Nor + e_Nar)),
    MII = log((Cyd + e_Cyd) / (Cox + e_Cox)),
    SRPI = (log(M00596_wt + 0.01) + log(Sred + e_Sred)) / 2,
    tea_mode = "ko"
  )
}

compute_indices_module_fallback <- function(comp_sub, w_vec) {
  wm <- function(mod) sum(get_col(comp_sub, mod) * w_vec, na.rm = TRUE)
  meth <- wm("M00567") + wm("M00357") + wm("M00356") + wm("M00563")
  ox <- wm("M00174")
  no3a <- wm("M00529")
  no3b <- wm("M00530")
  aer <- wm("M00155") + wm("M00154") + wm("M00416") + wm("M00417")
  mox <- wm("M00156") + wm("M00153")
  so4 <- wm("M00596")

  e <- function(x) max(x, 1e-6)
  data.table(
    MGI = log(meth + e(meth)),
    MFI = log((ox + 1e-6) / (meth + 1e-6)),
    DCI = log((no3b + 1e-6) / (no3a + 1e-6)),
    MII = log((mox + 1e-6) / (aer + 1e-6)),
    SRPI = log(so4 + 0.01),
    tea_mode = "module_completeness"
  )
}

sat_diag <- function(x, boundary_vals = BOUNDARY_VALS, eps = EPS) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(data.table(
      n = 0L,
      n_unique = 0L,
      sd = NA_real_,
      iqr = NA_real_,
      most_common_frac = NA_real_,
      max_boundary_frac = NA_real_,
      sat_flag = TRUE,
      sat_reason = "no_finite_values"
    ))
  }
  tx <- table(round(x, 6))
  most_common_frac <- max(tx) / n
  boundary_fracs <- sapply(boundary_vals, function(b) mean(abs(x - b) <= eps))
  max_boundary_frac <- max(boundary_fracs)
  n_unique <- uniqueN(round(x, 6))
  sdx <- sd(x)
  iqr_x <- IQR(x)

  sat_flag <- (n_unique <= 5L) || (!is.na(most_common_frac) && most_common_frac >= 0.8) ||
    (!is.na(max_boundary_frac) && max_boundary_frac >= 0.8) || (!is.na(sdx) && sdx < 1e-6)

  sat_reason <- if (n_unique <= 5L) {
    "few_unique_values"
  } else if (most_common_frac >= 0.8) {
    "single_value_dominant"
  } else if (max_boundary_frac >= 0.8) {
    "boundary_value_dominant"
  } else if (sdx < 1e-6) {
    "near_zero_sd"
  } else {
    "not_saturated"
  }

  data.table(
    n = as.integer(n),
    n_unique = as.integer(n_unique),
    sd = sdx,
    iqr = iqr_x,
    most_common_frac = most_common_frac,
    max_boundary_frac = max_boundary_frac,
    sat_flag = sat_flag,
    sat_reason = sat_reason
  )
}

compute_row_tea <- function(comp_wide, ko_wide, oap_per_taxon, taxa_vec, w_vec) {
  if (length(taxa_vec) == 0 || length(w_vec) == 0 || sum(w_vec, na.rm = TRUE) <= 0) {
    return(data.table(
      OAP = NA_real_, MGI = NA_real_, MFI = NA_real_, DCI = NA_real_, MII = NA_real_, SRPI = NA_real_,
      tea_mode = NA_character_,
      n_taxa_input = 0L,
      n_taxa_with_comp = 0L,
      n_taxa_with_ko = 0L,
      n_taxa_with_oap = 0L,
      weight_supported_comp = NA_real_,
      weight_supported_ko = NA_real_,
      weight_supported_oap = NA_real_
    ))
  }

  keep <- is.finite(w_vec) & (w_vec > 0)
  taxa_vec <- taxa_vec[keep]
  w_vec <- w_vec[keep]
  w_vec <- w_vec / sum(w_vec)

  comp_sub <- comp_wide[taxon %in% taxa_vec]
  if (nrow(comp_sub) == 0) {
    return(data.table(
      OAP = NA_real_, MGI = NA_real_, MFI = NA_real_, DCI = NA_real_, MII = NA_real_, SRPI = NA_real_,
      tea_mode = "no_annotation",
      n_taxa_input = length(taxa_vec),
      n_taxa_with_comp = 0L,
      n_taxa_with_ko = 0L,
      n_taxa_with_oap = 0L,
      weight_supported_comp = 0,
      weight_supported_ko = 0,
      weight_supported_oap = 0
    ))
  }

  setorder(comp_sub, taxon)
  ww <- w_vec[match(comp_sub$taxon, taxa_vec)]

  oap_join <- oap_per_taxon[taxon %in% comp_sub$taxon]
  OAP <- NA_real_
  if (nrow(oap_join) > 0) {
    oap_join[, w := ww[match(taxon, comp_sub$taxon)]]
    OAP <- weighted.mean(oap_join$oap_v3, oap_join$w, na.rm = TRUE)
  }

  ko_join <- ko_wide[taxon %in% comp_sub$taxon]
  if (nrow(ko_join) > 0) {
    ko_join[, w := ww[match(taxon, comp_sub$taxon)]]
    idx <- compute_indices_ko(ko_join, comp_sub, ww)
  } else {
    idx <- compute_indices_module_fallback(comp_sub, ww)
  }

  comp_taxa <- comp_sub$taxon
  ko_taxa <- if (nrow(ko_join) > 0) ko_join$taxon else character(0)
  oap_taxa <- if (nrow(oap_join) > 0) oap_join$taxon else character(0)

  data.table(
    OAP = OAP,
    MGI = idx$MGI,
    MFI = idx$MFI,
    DCI = idx$DCI,
    MII = idx$MII,
    SRPI = idx$SRPI,
    tea_mode = idx$tea_mode,
    n_taxa_input = length(taxa_vec),
    n_taxa_with_comp = length(comp_taxa),
    n_taxa_with_ko = uniqueN(ko_taxa),
    n_taxa_with_oap = uniqueN(oap_taxa),
    weight_supported_comp = sum(w_vec[taxa_vec %in% comp_taxa], na.rm = TRUE),
    weight_supported_ko = sum(w_vec[taxa_vec %in% ko_taxa], na.rm = TRUE),
    weight_supported_oap = sum(w_vec[taxa_vec %in% oap_taxa], na.rm = TRUE)
  )
}

center_within_sample <- function(module_long, exclude_grey = FALSE) {
  out <- copy(module_long)
  tag <- if (exclude_grey) "exclude_grey" else "include_grey"
  if (exclude_grey) out <- out[module != "grey"]
  out <- out[is.finite(value)]
  out[, sample_mean := mean(value, na.rm = TRUE), by = .(sample, index)]
  out[, centered := value - sample_mean]
  out[, analysis_set := tag]
  out
}

test_centered_vs_zero <- function(centered_dt, min_n = 20L, state_specific = FALSE) {
  grp_cols <- if (state_specific) c("state", "module", "index", "analysis_set") else c("module", "index", "analysis_set")
  out <- centered_dt[, {
    x <- centered[is.finite(centered)]
    n <- length(x)
    if (n < min_n) {
      .(n = n, median_centered = median(x, na.rm = TRUE), mean_centered = mean(x, na.rm = TRUE), p_value = NA_real_, status = "insufficient_n")
    } else {
      wt <- tryCatch(wilcox.test(x, mu = 0, exact = FALSE, conf.int = FALSE), error = function(e) e)
      if (inherits(wt, "error")) {
        .(n = n, median_centered = median(x, na.rm = TRUE), mean_centered = mean(x, na.rm = TRUE), p_value = NA_real_, status = "wilcox_failed")
      } else {
        .(n = n, median_centered = median(x, na.rm = TRUE), mean_centered = mean(x, na.rm = TRUE), p_value = wt$p.value, status = "ok")
      }
    }
  }, by = grp_cols]

  out[, q_value := p.adjust(p_value, method = "BH"), by = .(analysis_set, index)]
  out[, direction := fifelse(
    status != "ok" | !is.finite(q_value) | q_value >= 0.05,
    "no_call",
    fifelse(median_centered > 0, "enriched", fifelse(median_centered < 0, "depleted", "no_call"))
  )]
  out
}

log_msg("Loading main-window TEA inputs...")

meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
retained_taxa <- fread(file.path(DIRS$main, "taxa_after_filter.tsv"))
mods <- fread(file.path(DIRS$main, "module_assignments.tsv"))
hmm <- fread(file.path(DIRS$main, "hmm_states_main.tsv"))
func <- fread(INPUTS$prokaryote_function, select = c("taxon", "genome_id", "module_completeness_detail"))

tax <- if (grepl("\\.gz$", INPUTS$tax_damage)) {
  fread(cmd = sprintf("gzip -dc %s", shQuote(INPUTS$tax_damage)))
} else {
  fread(INPUTS$tax_damage)
}

sample_ids <- unique(meta_main$label)
retained <- unique(retained_taxa$taxon)

damaged_counts <- tax[
  is_dmg == "Damaged" &
    domain %in% c("d__Archaea", "d__Bacteria", "d__Viruses") &
    label %in% sample_ids &
    subspecies %in% retained,
  .(n_reads = sum(n_reads)),
  by = .(sample = label, taxon = subspecies)
]

sample_totals <- damaged_counts[, .(total_reads_retained_taxa = sum(n_reads)), by = sample]
damaged_counts <- merge(damaged_counts, sample_totals, by = "sample", all.x = TRUE)
damaged_counts[, rel_ab := n_reads / total_reads_retained_taxa]

if (any(!is.finite(damaged_counts$rel_ab))) {
  stop("Non-finite relative abundance detected in damaged_counts.")
}

missing_samples <- setdiff(sample_ids, unique(damaged_counts$sample))
if (length(missing_samples) > 0) {
  stop("Samples in main metadata with no retained damaged-read counts: ", paste(missing_samples, collapse = ", "))
}

mods <- mods[taxon %in% retained]
mods[is.na(module) | module == "", module := "unassigned"]

OAP_MODS <- rbindlist(lapply(names(TEA$oap_modules), function(m) {
  x <- TEA$oap_modules[[m]]
  data.table(module = m, tea_class = x$class, eo_mv = x$eo_mv, disc = x$disc)
}))

all_mods_needed <- unique(c(
  OAP_MODS$module,
  "M00174",
  "M00596",
  "M00567", "M00357", "M00356", "M00563",
  "M00529", "M00530",
  "M00155", "M00156", "M00153", "M00154", "M00416", "M00417", "M00973"
))

log_msg("Loading KEGG module completeness table for TEA annotations...")
kegg <- fread(
  cmd = sprintf("gzip -dc %s", shQuote(INPUTS$kegg_mods)),
  select = c("module", "genome_name", "stepwise_module_completeness", "enzyme_hits_in_module"),
  showProgress = FALSE
)
setnames(kegg, c("module", "reference", "completeness", "enzyme_hits"))
kegg <- kegg[module %in% all_mods_needed]

func[, ref_clean := genome_id]
func[grepl("^GCA_", genome_id), gcf_id := sub("^GCA_", "GCF_", genome_id)]

ref_lookup <- unique(rbind(
  func[, .(reference = ref_clean, taxon)],
  func[!is.na(gcf_id), .(reference = gcf_id, taxon)],
  fill = TRUE
)[!is.na(reference)])

kegg_taxa <- merge(kegg, ref_lookup, by = "reference", all.x = TRUE)
matched <- kegg_taxa[!is.na(taxon)]

fallback_taxa <- setdiff(func$taxon, matched$taxon)
if (length(fallback_taxa) > 0) {
  fb <- parse_modules_dt(func[taxon %in% fallback_taxa])
  fb <- fb[module %in% all_mods_needed]
  fb[, `:=`(reference = taxon, enzyme_hits = NA_character_)]
  taxon_mods <- rbind(
    matched[, .(module, taxon, completeness, enzyme_hits)],
    fb[, .(module, taxon, completeness, enzyme_hits)],
    fill = TRUE
  )
} else {
  taxon_mods <- matched[, .(module, taxon, completeness, enzyme_hits)]
}

taxon_mods <- taxon_mods[taxon %in% retained]
taxon_mods <- taxon_mods[is.finite(completeness)]

oap_taxon <- merge(
  taxon_mods[module %in% OAP_MODS$module],
  OAP_MODS,
  by = "module"
)[completeness >= PARAMS$oap_completeness_threshold,
  .(eff = max(completeness * disc, na.rm = TRUE)),
  by = .(taxon, tea_class, eo_mv)
]

oap_per_taxon <- oap_taxon[, .(
  oap_v3 = sum(eff * eo_mv) / sum(eff),
  n_classes = uniqueN(tea_class),
  dominant_class = tea_class[which.max(eff * eo_mv)],
  o2_eff = sum(eff[tea_class == "O2"], na.rm = TRUE),
  so4_eff = sum(eff[tea_class == "SO4"], na.rm = TRUE),
  ch4_eff = sum(eff[tea_class == "CH4"], na.rm = TRUE),
  no3_eff = sum(eff[tea_class == "NO3"], na.rm = TRUE)
), by = taxon]
oap_per_taxon[, oxy_class := fcase(
  o2_eff > 0 & so4_eff == 0 & ch4_eff == 0, "aerobic",
  so4_eff > 0 & o2_eff == 0, "sulfate_reducer",
  ch4_eff > 0 & o2_eff == 0, "methanogen",
  no3_eff > 0 & o2_eff == 0, "denitrifier",
  o2_eff > 0 & (so4_eff > 0 | ch4_eff > 0), "facultative",
  default = "no_call"
)]

ko_from_hits <- taxon_mods[!is.na(enzyme_hits) & enzyme_hits != "",
  .(taxon, module, enzyme_hits)
][, .(ko = unlist(strsplit(enzyme_hits, ",", fixed = TRUE))), by = .(taxon, module)
][ko %in% TEA$target_kos]

ko_wide <- if (nrow(ko_from_hits) > 0) {
  dcast(ko_from_hits[, .N, by = .(taxon, ko)], taxon ~ ko, value.var = "N", fill = 0)
} else {
  data.table(taxon = character())
}

best_comp <- taxon_mods[, .(best_comp = max(completeness, na.rm = TRUE)), by = .(taxon, module)]
comp_wide <- dcast(best_comp, taxon ~ module, value.var = "best_comp", fill = 0)

log_msg("Computing sample-level TEA on main-window set...")
sample_rows <- lapply(sample_ids, function(samp) {
  sub <- damaged_counts[sample == samp]
  out <- compute_row_tea(
    comp_wide = comp_wide,
    ko_wide = ko_wide,
    oap_per_taxon = oap_per_taxon,
    taxa_vec = sub$taxon,
    w_vec = sub$rel_ab
  )
  out[, sample := samp]
  out[, total_reads_retained_taxa := unique(sub$total_reads_retained_taxa)]
  out
})
tea_sample <- rbindlist(sample_rows, fill = TRUE)
setcolorder(tea_sample, c("sample", "total_reads_retained_taxa", INDEX_COLS, setdiff(names(tea_sample), c("sample", "total_reads_retained_taxa", INDEX_COLS))))

log_msg("Computing sample x module TEA (grey retained)...")
counts_mod <- merge(damaged_counts, mods, by = "taxon", all.x = TRUE)
counts_mod[is.na(module) | module == "", module := "unassigned"]

module_rows <- counts_mod[, {
  module_abundance <- sum(rel_ab, na.rm = TRUE)
  w_mod <- rel_ab / module_abundance
  rr <- compute_row_tea(
    comp_wide = comp_wide,
    ko_wide = ko_wide,
    oap_per_taxon = oap_per_taxon,
    taxa_vec = taxon,
    w_vec = w_mod
  )
  rr[, `:=`(
    module_abundance = module_abundance,
    n_taxa_in_module_sample = uniqueN(taxon),
    module_total_reads = sum(n_reads)
  )]
  rr
}, by = .(sample, module)]

tea_sample_module <- module_rows[module_abundance > 0]
setcolorder(tea_sample_module, c("sample", "module", "module_abundance", INDEX_COLS, setdiff(names(tea_sample_module), c("sample", "module", "module_abundance", INDEX_COLS))))

hmm_use <- unique(hmm[, .(sample, core, age_kyr, state, state_label = label)])
meta_use <- unique(meta_main[, .(sample = label, core, age_kyr)])

tea_sample <- merge(tea_sample, meta_use, by = "sample", all.x = TRUE)
tea_sample <- merge(tea_sample, hmm_use[, .(sample, state, state_label)], by = "sample", all.x = TRUE)

tea_sample_module <- merge(tea_sample_module, meta_use, by = "sample", all.x = TRUE)
tea_sample_module <- merge(tea_sample_module, hmm_use[, .(sample, state, state_label)], by = "sample", all.x = TRUE)

log_msg("Building TEA summaries and conservative tests...")

module_summary <- tea_sample_module[, .(
  n_sample_rows = .N,
  n_samples = uniqueN(sample),
  mean_module_abundance = mean(module_abundance, na.rm = TRUE),
  median_module_abundance = median(module_abundance, na.rm = TRUE),
  frac_ko_mode = mean(tea_mode == "ko", na.rm = TRUE),
  mean_weight_supported_comp = mean(weight_supported_comp, na.rm = TRUE),
  mean_weight_supported_ko = mean(weight_supported_ko, na.rm = TRUE),
  mean_weight_supported_oap = mean(weight_supported_oap, na.rm = TRUE),
  mean_OAP = mean(OAP, na.rm = TRUE),
  median_OAP = median(OAP, na.rm = TRUE),
  mean_MGI = mean(MGI, na.rm = TRUE),
  median_MGI = median(MGI, na.rm = TRUE),
  mean_MFI = mean(MFI, na.rm = TRUE),
  median_MFI = median(MFI, na.rm = TRUE),
  mean_DCI = mean(DCI, na.rm = TRUE),
  median_DCI = median(DCI, na.rm = TRUE),
  mean_MII = mean(MII, na.rm = TRUE),
  median_MII = median(MII, na.rm = TRUE),
  mean_SRPI = mean(SRPI, na.rm = TRUE),
  median_SRPI = median(SRPI, na.rm = TRUE)
), by = module][order(module)]

state_summary <- melt(
  tea_sample,
  id.vars = c("sample", "core", "age_kyr", "state", "state_label", "tea_mode"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)[is.finite(value), .(
  n_samples = .N,
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  q25 = quantile(value, 0.25, na.rm = TRUE),
  q75 = quantile(value, 0.75, na.rm = TRUE)
), by = .(state, state_label, index)][order(index, state)]

sample_long <- melt(
  tea_sample,
  id.vars = c("sample", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)[is.finite(value)]

sample_state_tests <- sample_long[, {
  s_counts <- .SD[, .N, by = state]
  enough <- nrow(s_counts) >= 2 && all(s_counts$N >= MIN_N_STATE_EFFECT_PER_STATE)
  if (!enough) {
    .(n = .N, n_states = uniqueN(state), kruskal_p = NA_real_, lm_coreadj_p = NA_real_, status = "insufficient_by_state")
  } else {
    kt <- tryCatch(kruskal.test(value ~ factor(state), data = .SD), error = function(e) e)
    lmfit <- tryCatch(lm(value ~ factor(state) + factor(core), data = .SD), error = function(e) e)
    p_lm <- NA_real_
    if (!inherits(lmfit, "error")) {
      an <- anova(lmfit)
      rn <- rownames(an)
      idx_term <- grep("factor\\(state\\)", rn)
      if (length(idx_term) == 1) p_lm <- an[idx_term, "Pr(>F)"]
    }
    .(
      n = .N,
      n_states = uniqueN(state),
      kruskal_p = if (!inherits(kt, "error")) kt$p.value else NA_real_,
      lm_coreadj_p = p_lm,
      status = "ok"
    )
  }
}, by = index]
sample_state_tests[, kruskal_q := p.adjust(kruskal_p, method = "BH")]
sample_state_tests[, lm_coreadj_q := p.adjust(lm_coreadj_p, method = "BH")]

module_long <- melt(
  tea_sample_module,
  id.vars = c("sample", "module", "module_abundance", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)

center_inc <- center_within_sample(module_long, exclude_grey = FALSE)
center_exc <- center_within_sample(module_long, exclude_grey = TRUE)
center_all <- rbindlist(list(center_inc, center_exc), fill = TRUE)

module_enrichment <- test_centered_vs_zero(center_all, min_n = MIN_N_MODULE_TEST, state_specific = FALSE)

state_module_enrichment <- test_centered_vs_zero(
  center_all[!is.na(state)],
  min_n = MIN_N_STATE_MODULE_TEST,
  state_specific = TRUE
)

state_effect <- center_inc[!is.na(state), {
  s_counts <- .SD[, .N, by = state]
  enough <- nrow(s_counts) >= 2 && all(s_counts$N >= MIN_N_STATE_EFFECT_PER_STATE)
  if (!enough) {
    .(n = .N, n_states = uniqueN(state), kruskal_p = NA_real_, lm_coreadj_p = NA_real_, status = "insufficient_by_state")
  } else {
    kt <- tryCatch(kruskal.test(centered ~ factor(state), data = .SD), error = function(e) e)
    lmfit <- tryCatch(lm(centered ~ factor(state) + factor(core), data = .SD), error = function(e) e)
    p_lm <- NA_real_
    if (!inherits(lmfit, "error")) {
      an <- anova(lmfit)
      rn <- rownames(an)
      idx_term <- grep("factor\\(state\\)", rn)
      if (length(idx_term) == 1) p_lm <- an[idx_term, "Pr(>F)"]
    }
    .(
      n = .N,
      n_states = uniqueN(state),
      kruskal_p = if (!inherits(kt, "error")) kt$p.value else NA_real_,
      lm_coreadj_p = p_lm,
      status = "ok"
    )
  }
}, by = .(module, index)]
state_effect[, kruskal_q := p.adjust(kruskal_p, method = "BH")]
state_effect[, lm_coreadj_q := p.adjust(lm_coreadj_p, method = "BH")]

sat_sample <- rbindlist(lapply(INDEX_COLS, function(idx) {
  s <- sat_diag(tea_sample[[idx]])
  s[, `:=`(scope = "sample", group = "all", index = idx)]
  s
}), fill = TRUE)

sat_module <- rbindlist(lapply(INDEX_COLS, function(idx) {
  tea_sample_module[, {
    s <- sat_diag(get(idx))
    s[, `:=`(scope = "sample_module", group = as.character(module), index = idx)]
    s
  }, by = module]
}), fill = TRUE)

sat_centered <- rbindlist(lapply(INDEX_COLS, function(idx) {
  center_inc[index == idx, {
    s <- sat_diag(centered)
    s[, `:=`(scope = "centered_include_grey", group = as.character(module), index = idx)]
    s
  }, by = module]
}), fill = TRUE)

saturation <- rbindlist(list(sat_sample, sat_module, sat_centered), fill = TRUE)

taxon_annot <- unique(data.table(taxon = retained))
taxon_annot <- merge(taxon_annot, mods[, .(taxon, module)], by = "taxon", all.x = TRUE)
taxon_annot <- merge(
  taxon_annot,
  comp_wide[, .(taxon, has_comp = TRUE)],
  by = "taxon",
  all.x = TRUE
)
taxon_annot <- merge(
  taxon_annot,
  ko_wide[, .(taxon, has_ko = TRUE)],
  by = "taxon",
  all.x = TRUE
)
taxon_annot <- merge(
  taxon_annot,
  oap_per_taxon[, .(taxon, oap_v3, n_classes, dominant_class, oxy_class, o2_eff, so4_eff, ch4_eff, no3_eff, has_oap = TRUE)],
  by = "taxon",
  all.x = TRUE
)
for (v in c("has_comp", "has_ko", "has_oap")) taxon_annot[is.na(get(v)), (v) := FALSE]

ko_presence_long <- ko_from_hits[, .(taxon, ko)][, .(has_ko = TRUE), by = .(taxon, ko)]
ko_presence <- if (nrow(ko_presence_long) > 0) {
  dcast(ko_presence_long, taxon ~ ko, value.var = "has_ko", fill = FALSE)
} else {
  data.table(taxon = retained)
}
for (ko in setdiff(TEA$target_kos, names(ko_presence))) ko_presence[, (ko) := FALSE]
for (ko in TEA$target_kos) ko_presence[, (ko) := as.logical(get(ko))]

ko_flags <- ko_presence[, .(
  taxon,
  has_mcr = K00399 | K00401 | K00402,
  has_pmo = K14080,
  has_nar = K00370,
  has_nos = K00376,
  has_cyd = K00425,
  has_cox = K02274 | K02276,
  has_sred = K00394 | K00395 | K00958 | K11180 | K11181
)]
taxon_annot <- merge(taxon_annot, ko_flags, by = "taxon", all.x = TRUE)
for (v in c("has_mcr", "has_pmo", "has_nar", "has_nos", "has_cyd", "has_cox", "has_sred")) {
  taxon_annot[is.na(get(v)), (v) := FALSE]
}

coverage_sample <- tea_sample[, .(
  n_taxa_input,
  n_taxa_with_comp,
  n_taxa_with_ko,
  n_taxa_with_oap,
  frac_taxa_with_comp = n_taxa_with_comp / pmax(n_taxa_input, 1),
  frac_taxa_with_ko = n_taxa_with_ko / pmax(n_taxa_input, 1),
  frac_taxa_with_oap = n_taxa_with_oap / pmax(n_taxa_input, 1),
  weight_supported_comp,
  weight_supported_ko,
  weight_supported_oap,
  tea_mode,
  state,
  state_label,
  core,
  age_kyr
), by = sample]

coverage_module <- tea_sample_module[, .(
  sample,
  module,
  module_abundance,
  n_taxa_in_module_sample,
  n_taxa_with_comp,
  n_taxa_with_ko,
  n_taxa_with_oap,
  frac_taxa_with_comp = n_taxa_with_comp / pmax(n_taxa_in_module_sample, 1),
  frac_taxa_with_ko = n_taxa_with_ko / pmax(n_taxa_in_module_sample, 1),
  frac_taxa_with_oap = n_taxa_with_oap / pmax(n_taxa_in_module_sample, 1),
  weight_supported_comp,
  weight_supported_ko,
  weight_supported_oap,
  tea_mode,
  state,
  state_label,
  core,
  age_kyr
)]

coverage_summary <- data.table(
  metric = c(
    "n_samples_main",
    "n_retained_taxa",
    "n_taxa_with_module_assignment",
    "n_taxa_with_any_comp",
    "n_taxa_with_any_ko",
    "n_taxa_with_any_oap",
    "median_sample_weight_supported_comp",
    "median_sample_weight_supported_ko",
    "median_sample_weight_supported_oap",
    "median_module_weight_supported_comp",
    "median_module_weight_supported_ko",
    "median_module_weight_supported_oap"
  ),
  value = c(
    nrow(tea_sample),
    length(retained),
    nrow(taxon_annot[!is.na(module)]),
    nrow(taxon_annot[has_comp == TRUE]),
    nrow(taxon_annot[has_ko == TRUE]),
    nrow(taxon_annot[has_oap == TRUE]),
    median(coverage_sample$weight_supported_comp, na.rm = TRUE),
    median(coverage_sample$weight_supported_ko, na.rm = TRUE),
    median(coverage_sample$weight_supported_oap, na.rm = TRUE),
    median(coverage_module$weight_supported_comp, na.rm = TRUE),
    median(coverage_module$weight_supported_ko, na.rm = TRUE),
    median(coverage_module$weight_supported_oap, na.rm = TRUE)
  )
)

decision_notes <- rbindlist(list(
  data.table(note_key = "scope", note_value = "Main-window only (<=150 ka), using existing clean workflow outputs for retained taxa/modules/HMM states."),
  data.table(note_key = "abundance_weighting", note_value = "Damaged-read relative abundance computed across retained taxa per sample; module-level TEA uses within-sample, within-module normalized weights."),
  data.table(note_key = "tea_mode", note_value = "KO-first TEA logic retained; module-completeness fallback used when KO evidence unavailable."),
  data.table(note_key = "grey_policy", note_value = "Primary module centering and enrichment tests include grey; exclude-grey sensitivity outputs are also exported."),
  data.table(note_key = "interpretation", note_value = "TEA indices represent annotation-weighted genomic potential, not direct in situ activity."),
  data.table(note_key = "state_reuse", note_value = "HMM states are reused from 03_hmm_main outputs; states are not refit in this TEA stage."),
  data.table(note_key = "taxon_trait_scope", note_value = "Taxon-level OAP and metabolic calls are annotation-derived trait proxies (genomic potential), not measured metabolic activity.")
), fill = TRUE)

fwrite(tea_sample, file.path(DIRS$main_tea, "tea_indices_per_sample_main.tsv"), sep = "\t")
fwrite(tea_sample_module, file.path(DIRS$main_tea, "tea_indices_per_sample_module_main.tsv"), sep = "\t")
fwrite(module_summary, file.path(DIRS$main_tea, "tea_module_summary_main.tsv"), sep = "\t")
fwrite(state_summary, file.path(DIRS$main_tea, "tea_state_summary_main.tsv"), sep = "\t")
fwrite(sample_state_tests, file.path(DIRS$main_tea, "tea_sample_level_state_tests.tsv"), sep = "\t")
fwrite(module_enrichment, file.path(DIRS$main_tea, "tea_module_enrichment_centered_wilcoxon.tsv"), sep = "\t")
fwrite(state_module_enrichment, file.path(DIRS$main_tea, "tea_state_module_enrichment_centered_wilcoxon.tsv"), sep = "\t")
fwrite(state_effect, file.path(DIRS$main_tea, "tea_module_index_state_effects.tsv"), sep = "\t")
fwrite(center_all, file.path(DIRS$main_tea, "tea_module_centered_values.tsv"), sep = "\t")
fwrite(saturation, file.path(DIRS$main_tea, "tea_saturation_diagnostics.tsv"), sep = "\t")
fwrite(coverage_sample, file.path(DIRS$main_tea, "tea_annotation_coverage_sample.tsv"), sep = "\t")
fwrite(coverage_module, file.path(DIRS$main_tea, "tea_annotation_coverage_sample_module.tsv"), sep = "\t")
fwrite(taxon_annot, file.path(DIRS$main_tea, "tea_taxon_metabolic_calls.tsv"), sep = "\t")
fwrite(taxon_annot, file.path(DIRS$main_tea, "tea_annotation_coverage_taxon.tsv"), sep = "\t")
fwrite(coverage_summary, file.path(DIRS$main_tea, "tea_annotation_coverage_summary.tsv"), sep = "\t")
fwrite(decision_notes, file.path(DIRS$main_tea, "tea_decision_notes.tsv"), sep = "\t")

write_run_metadata(
  file.path(DIRS$main_tea, "04a_tea_main_run_metadata.tsv"),
  "04a_tea_main.R",
  extra = list(
    sample_scope = "main_window_training_plus_validation",
    max_age_kyr = PARAMS$main_max_age_kyr,
    n_samples = nrow(tea_sample),
    n_sample_module_rows = nrow(tea_sample_module),
    n_taxa_annotation_rows = nrow(taxon_annot),
    n_taxa_with_oap = nrow(taxon_annot[has_oap == TRUE]),
    oap_threshold = PARAMS$oap_completeness_threshold,
    module_test_min_n = MIN_N_MODULE_TEST,
    state_module_test_min_n = MIN_N_STATE_MODULE_TEST,
    state_effect_min_n_per_state = MIN_N_STATE_EFFECT_PER_STATE,
    include_grey_primary = TRUE,
    exclude_grey_sensitivity = TRUE,
    interpretation_scope = "genomic_potential_not_direct_activity"
  )
)
write_session_info(file.path(DIRS$main_tea, "04a_tea_main_sessionInfo.txt"))

log_msg("04a_tea_main complete.")
