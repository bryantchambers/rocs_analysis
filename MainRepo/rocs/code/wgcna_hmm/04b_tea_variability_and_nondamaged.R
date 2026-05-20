#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)

OUT_DIR <- file.path(DIRS$main, "tea_variability")
ensure_dirs(c(OUT_DIR))

INDEX_COLS <- c("OAP", "MGI", "MFI", "DCI", "MII", "SRPI")

required_inputs <- c(
  INPUTS$tax_damage,
  INPUTS$kegg_mods,
  INPUTS$prokaryote_function,
  file.path(DIRS$main, "sample_metadata_main.tsv"),
  file.path(DIRS$main, "module_assignments.tsv"),
  file.path(DIRS$main, "hmm_states_main.tsv"),
  file.path(DIRS$main_tea, "tea_indices_per_sample_main.tsv"),
  file.path(DIRS$main_tea, "tea_indices_per_sample_module_main.tsv"),
  file.path(DIRS$main_tea, "tea_module_centered_values.tsv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required 04b inputs: ", paste(missing_inputs, collapse = "; "))
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

cv_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  mu <- mean(x)
  if (!is.finite(mu) || abs(mu) < 1e-8) return(NA_real_)
  stats::sd(x) / abs(mu)
}

iqr_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  IQR(x)
}

robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  madv <- mad(x, center = med, constant = 1, na.rm = TRUE)
  if (!is.finite(madv) || madv < 1e-8) {
    return(rep(NA_real_, length(x)))
  }
  (x - med) / madv
}

pair_metrics <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) == 0) {
    return(data.table(
      n_indices = 0L,
      euclidean = NA_real_,
      manhattan = NA_real_,
      cosine_similarity = NA_real_,
      pearson = NA_real_,
      spearman = NA_real_
    ))
  }
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  cos_sim <- if (denom > 0) sum(x * y) / denom else NA_real_
  pcor <- if (length(unique(x)) > 1 && length(unique(y)) > 1) suppressWarnings(cor(x, y, method = "pearson")) else NA_real_
  scor <- if (length(unique(x)) > 1 && length(unique(y)) > 1) suppressWarnings(cor(x, y, method = "spearman")) else NA_real_
  data.table(
    n_indices = as.integer(length(x)),
    euclidean = sqrt(sum((x - y)^2)),
    manhattan = sum(abs(x - y)),
    cosine_similarity = cos_sim,
    pearson = pcor,
    spearman = scor
  )
}

log_msg("Loading TEA outputs for variability/specialization summaries...")

meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
hmm <- fread(file.path(DIRS$main, "hmm_states_main.tsv"))
tea_sample <- fread(file.path(DIRS$main_tea, "tea_indices_per_sample_main.tsv"))
tea_sample_module <- fread(file.path(DIRS$main_tea, "tea_indices_per_sample_module_main.tsv"))

module_long <- melt(
  tea_sample_module,
  id.vars = c("sample", "module", "module_abundance", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)
module_long <- module_long[is.finite(value)]
module_long[, sample_index_mean := mean(value, na.rm = TRUE), by = .(sample, index)]
module_long[, centered_within_sample := value - sample_index_mean]

module_variability <- module_long[, .(
  n_rows = .N,
  n_samples = uniqueN(sample),
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  iqr = iqr_safe(value),
  cv = cv_safe(value),
  centered_mean = mean(centered_within_sample, na.rm = TRUE),
  centered_median = median(centered_within_sample, na.rm = TRUE),
  centered_abs_mean = mean(abs(centered_within_sample), na.rm = TRUE),
  mean_module_abundance = mean(module_abundance, na.rm = TRUE),
  median_module_abundance = median(module_abundance, na.rm = TRUE)
), by = .(module, index)][order(index, module)]

state_module_variability <- module_long[, .(
  n_rows = .N,
  n_samples = uniqueN(sample),
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  iqr = iqr_safe(value),
  cv = cv_safe(value),
  centered_mean = mean(centered_within_sample, na.rm = TRUE),
  centered_median = median(centered_within_sample, na.rm = TRUE),
  specialization_abs_centered_mean = mean(abs(centered_within_sample), na.rm = TRUE),
  mean_module_abundance = mean(module_abundance, na.rm = TRUE)
), by = .(state, state_label, module, index)][order(index, state, module)]

log_msg("Computing conservative temporal descriptive summaries...")

temporal_stats <- function(dt, value_col = "value") {
  x <- dt[[value_col]]
  y <- dt$age_kyr
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  rho <- NA_real_
  p <- NA_real_
  if (n >= 6 && length(unique(y)) >= 4 && length(unique(x)) >= 3) {
    ct <- tryCatch(cor.test(y, x, method = "spearman", exact = FALSE), error = function(e) e)
    if (!inherits(ct, "error")) {
      rho <- unname(ct$estimate)
      p <- ct$p.value
    }
  }
  data.table(
    n = n,
    age_min_kyr = if (n > 0) min(y) else NA_real_,
    age_max_kyr = if (n > 0) max(y) else NA_real_,
    mean = if (n > 0) mean(x) else NA_real_,
    median = if (n > 0) median(x) else NA_real_,
    sd = if (n > 1) sd(x) else NA_real_,
    iqr = if (n > 0) IQR(x) else NA_real_,
    spearman_rho_age = rho,
    spearman_p_age = p
  )
}

temporal_module_core <- module_long[, temporal_stats(.SD), by = .(core, module, index)]
temporal_module_core[, `:=`(scope = "sample_module", group_level = "core")]

temporal_module_state <- module_long[, temporal_stats(.SD), by = .(state, state_label, module, index)]
temporal_module_state[, `:=`(scope = "sample_module", group_level = "state")]

sample_long <- melt(
  tea_sample,
  id.vars = c("sample", "core", "age_kyr", "state", "state_label"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)
sample_long <- sample_long[is.finite(value)]

temporal_sample_core <- sample_long[, temporal_stats(.SD), by = .(core, index)]
temporal_sample_core[, `:=`(scope = "sample", group_level = "core")]

temporal_sample_state <- sample_long[, temporal_stats(.SD), by = .(state, state_label, index)]
temporal_sample_state[, `:=`(scope = "sample", group_level = "state")]

temporal_summary <- rbindlist(
  list(temporal_module_core, temporal_module_state, temporal_sample_core, temporal_sample_state),
  fill = TRUE
)
temporal_summary[, spearman_q_age := p.adjust(spearman_p_age, method = "BH"), by = .(scope, group_level, index)]

log_msg("Calling dominant/suppressed indices (descriptive only)...")

module_calls <- copy(module_variability)
module_calls[, robust_z_centered_median := robust_z(centered_median), by = index]
module_calls[, dominance_call := fcase(
  is.finite(robust_z_centered_median) & robust_z_centered_median >= 0.75, "relatively_dominant",
  is.finite(robust_z_centered_median) & robust_z_centered_median <= -0.75, "relatively_suppressed",
  default = "intermediate"
)]
module_calls[, call_scope := "module"]

state_values <- sample_long[, .(
  n_samples = .N,
  state_mean = mean(value, na.rm = TRUE),
  state_median = median(value, na.rm = TRUE)
), by = .(state, state_label, index)]
state_values[, robust_z_state_median := robust_z(state_median), by = index]
state_values[, dominance_call := fcase(
  is.finite(robust_z_state_median) & robust_z_state_median >= 0.75, "relatively_high",
  is.finite(robust_z_state_median) & robust_z_state_median <= -0.75, "relatively_low",
  default = "intermediate"
)]
state_values[, call_scope := "state"]

dominance_calls <- rbindlist(list(
  module_calls[, .(call_scope, module, state = NA_integer_, state_label = NA_character_, index, n = n_rows, reference_median = centered_median, robust_z = robust_z_centered_median, dominance_call)],
  state_values[, .(call_scope, module = NA_character_, state, state_label, index, n = n_samples, reference_median = state_median, robust_z = robust_z_state_median, dominance_call)]
), fill = TRUE)

log_msg("Building non-damaged TEA profiles on the same main sample set...")

tax <- if (grepl("\\.gz$", INPUTS$tax_damage)) {
  fread(cmd = sprintf("gzip -dc %s", shQuote(INPUTS$tax_damage)))
} else {
  fread(INPUTS$tax_damage)
}

func <- fread(INPUTS$prokaryote_function, select = c("taxon", "genome_id", "module_completeness_detail"))

sample_ids <- unique(meta_main$label)

nondmg_counts <- tax[
  is_dmg == "Non-damaged" &
    label %in% sample_ids &
    !is.na(subspecies) & subspecies != "",
  .(n_reads = sum(n_reads)),
  by = .(sample = label, taxon = subspecies)
]

if (nrow(nondmg_counts) == 0) {
  stop("No non-damaged reads found in main sample set; cannot compute non-damaged TEA summaries.")
}

# Conservative default in this run: no additional filter.
nondmg_counts[, passed_light_filter := TRUE]

nondmg_sample_totals <- nondmg_counts[, .(
  nondmg_total_reads = sum(n_reads),
  nondmg_n_taxa = uniqueN(taxon)
), by = sample]
nondmg_counts <- merge(nondmg_counts, nondmg_sample_totals, by = "sample", all.x = TRUE)
nondmg_counts[, rel_ab := n_reads / nondmg_total_reads]

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

nondmg_taxa <- unique(nondmg_counts$taxon)
taxon_mods <- taxon_mods[taxon %in% nondmg_taxa & is.finite(completeness)]

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

nondmg_sample_rows <- lapply(sample_ids, function(samp) {
  sub <- nondmg_counts[sample == samp]
  rr <- compute_row_tea(
    comp_wide = comp_wide,
    ko_wide = ko_wide,
    oap_per_taxon = oap_per_taxon,
    taxa_vec = sub$taxon,
    w_vec = sub$rel_ab
  )
  rr[, sample := samp]
  rr[, nondmg_total_reads := if (nrow(sub) > 0) unique(sub$nondmg_total_reads) else 0]
  rr[, nondmg_n_taxa := if (nrow(sub) > 0) unique(sub$nondmg_n_taxa) else 0]
  rr
})
nondmg_sample <- rbindlist(nondmg_sample_rows, fill = TRUE)

meta_use <- unique(meta_main[, .(sample = label, core, age_kyr)])
hmm_use <- unique(hmm[, .(sample, state, state_label = label)])
nondmg_sample <- merge(nondmg_sample, meta_use, by = "sample", all.x = TRUE)
nondmg_sample <- merge(nondmg_sample, hmm_use, by = "sample", all.x = TRUE)

nondmg_sample_long <- melt(
  nondmg_sample,
  id.vars = c("sample", "core", "age_kyr", "state", "state_label", "nondmg_total_reads", "nondmg_n_taxa"),
  measure.vars = INDEX_COLS,
  variable.name = "index",
  value.name = "value"
)[is.finite(value)]

nondmg_aggregate <- nondmg_sample_long[, .(
  n_samples = .N,
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  iqr = IQR(value, na.rm = TRUE),
  read_weighted_mean = weighted.mean(value, w = nondmg_total_reads, na.rm = TRUE)
), by = index][order(index)]

nondmg_taxon_stats <- nondmg_counts[, .(
  n_samples_non_damaged = sum(n_reads > 0),
  total_reads_non_damaged = sum(n_reads)
), by = taxon]

nondmg_taxon_annot <- unique(data.table(taxon = nondmg_taxa))
nondmg_taxon_annot <- merge(nondmg_taxon_annot, nondmg_taxon_stats, by = "taxon", all.x = TRUE)
nondmg_taxon_annot <- merge(nondmg_taxon_annot, comp_wide[, .(taxon, has_comp = TRUE)], by = "taxon", all.x = TRUE)
nondmg_taxon_annot <- merge(nondmg_taxon_annot, ko_wide[, .(taxon, has_ko = TRUE)], by = "taxon", all.x = TRUE)
nondmg_taxon_annot <- merge(
  nondmg_taxon_annot,
  oap_per_taxon[, .(taxon, oap_v3, n_classes, dominant_class, o2_eff, so4_eff, ch4_eff, no3_eff, has_oap = TRUE)],
  by = "taxon",
  all.x = TRUE
)
for (v in c("has_comp", "has_ko", "has_oap")) nondmg_taxon_annot[is.na(get(v)), (v) := FALSE]
nondmg_taxon_annot[, light_filter_rule := "none (full Non-damaged taxa in main sample set retained)"]

log_msg("Comparing non-damaged aggregate TEA signatures to damaged modules...")

module_agg <- module_long[, .(
  module_mean = mean(value, na.rm = TRUE),
  module_median = median(value, na.rm = TRUE)
), by = .(module, index)]

nd_vec <- nondmg_aggregate$mean
names(nd_vec) <- nondmg_aggregate$index

module_comp <- rbindlist(lapply(unique(module_agg$module), function(m) {
  sub <- module_agg[module == m]
  mod_vec <- sub$module_mean
  names(mod_vec) <- sub$index
  idx <- intersect(names(nd_vec), names(mod_vec))
  mm <- pair_metrics(nd_vec[idx], mod_vec[idx])
  mm[, module := m]
  mm
}), fill = TRUE)

module_comp[, rank_by_euclidean := frank(euclidean, ties.method = "min")]
module_comp <- module_comp[order(rank_by_euclidean, module)]

nondmg_vs_module_by_index <- merge(
  nondmg_aggregate[, .(index, nondmg_mean = mean, nondmg_median = median)],
  module_agg,
  by = "index",
  all = TRUE
)
nondmg_vs_module_by_index[, `:=`(
  diff_mean = module_mean - nondmg_mean,
  abs_diff_mean = abs(module_mean - nondmg_mean)
)]

brown_idx <- nondmg_vs_module_by_index[module == "brown"]
brown_idx[, rank_abs_diff_among_modules := frank(abs_diff_mean, ties.method = "min"), by = index]

brown_sample <- tea_sample_module[module == "brown", c("sample", INDEX_COLS), with = FALSE]
setnames(brown_sample, INDEX_COLS, paste0("brown_", INDEX_COLS))

nd_sample_cmp <- merge(
  nondmg_sample[, c("sample", INDEX_COLS), with = FALSE],
  brown_sample,
  by = "sample",
  all = FALSE
)

brown_sample_match <- rbindlist(lapply(INDEX_COLS, function(idx) {
  x <- nd_sample_cmp[[idx]]
  y <- nd_sample_cmp[[paste0("brown_", idx)]]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  data.table(
    comparison_scope = "sample_matched",
    index = idx,
    n_samples = n,
    mean_diff_brown_minus_nondmg = if (n > 0) mean(y - x) else NA_real_,
    median_diff_brown_minus_nondmg = if (n > 0) median(y - x) else NA_real_,
    rmse = if (n > 0) sqrt(mean((y - x)^2)) else NA_real_,
    pearson = if (n >= 3 && length(unique(x)) > 1 && length(unique(y)) > 1) suppressWarnings(cor(y, x, method = "pearson")) else NA_real_,
    spearman = if (n >= 3 && length(unique(x)) > 1 && length(unique(y)) > 1) suppressWarnings(cor(y, x, method = "spearman")) else NA_real_
  )
}), fill = TRUE)

brown_global <- pair_metrics(
  brown_idx$nondmg_mean,
  brown_idx$module_mean
)
brown_global[, `:=`(
  comparison_scope = "aggregate_vector",
  module = "brown",
  rank_by_euclidean_among_modules = module_comp[module == "brown", rank_by_euclidean]
)]

brown_vs_nondmg <- rbindlist(list(
  brown_idx[, .(
    comparison_scope = "per_index_aggregate",
    index,
    nondmg_mean,
    nondmg_median,
    brown_mean = module_mean,
    brown_median = module_median,
    diff_mean_brown_minus_nondmg = diff_mean,
    abs_diff_mean,
    rank_abs_diff_among_modules
  )],
  brown_sample_match,
  brown_global[, .(
    comparison_scope,
    index = NA_character_,
    nondmg_mean = NA_real_,
    nondmg_median = NA_real_,
    brown_mean = NA_real_,
    brown_median = NA_real_,
    diff_mean_brown_minus_nondmg = NA_real_,
    abs_diff_mean = NA_real_,
    rank_abs_diff_among_modules = as.numeric(rank_by_euclidean_among_modules),
    n_indices,
    euclidean,
    manhattan,
    cosine_similarity,
    pearson,
    spearman
  )]
), fill = TRUE)

decision_notes <- rbindlist(list(
  data.table(note_key = "scope", note_value = "04b extends 04a with descriptive TEA variability/specialization, temporal summaries, and non-damaged profile comparisons in main-window samples (<=150 ka)."),
  data.table(note_key = "nondamaged_taxa_scope", note_value = "Used full set of taxa with is_dmg == Non-damaged in main samples; no reuse of damaged-only taxa_after_filter."),
  data.table(note_key = "nondamaged_filter", note_value = "No additional light filter applied in this run (conservative full-set handling)."),
  data.table(note_key = "transformations", note_value = "No CLR or additional compositional transformation used for non-damaged TEA comparisons."),
  data.table(note_key = "comparison_interpretation", note_value = "Non-damaged vs damaged-module comparisons are descriptive TEA phenotype/profile comparisons, not network membership inference."),
  data.table(note_key = "dominance_logic", note_value = "Dominance/suppression calls are relative within-index robust-z descriptive calls; not formal hypothesis tests of biological activity."),
  data.table(note_key = "temporal_logic", note_value = "Temporal analysis uses conservative descriptive summaries and Spearman age association by core/state; no complex time-series modeling."),
  data.table(note_key = "interpretation", note_value = "TEA indices are annotation-weighted genomic potential proxies, not direct measured in situ activity.")
), fill = TRUE)

fwrite(module_variability, file.path(OUT_DIR, "tea_module_variability_summary.tsv"), sep = "\t")
fwrite(state_module_variability, file.path(OUT_DIR, "tea_state_module_variability_specialization_summary.tsv"), sep = "\t")
fwrite(temporal_summary, file.path(OUT_DIR, "tea_temporal_descriptive_summary.tsv"), sep = "\t")
fwrite(dominance_calls, file.path(OUT_DIR, "tea_dominant_metabolism_calls.tsv"), sep = "\t")

fwrite(nondmg_taxon_annot, file.path(OUT_DIR, "tea_nondamaged_taxon_annotation.tsv"), sep = "\t")
fwrite(nondmg_sample, file.path(OUT_DIR, "tea_nondamaged_indices_per_sample_main.tsv"), sep = "\t")
fwrite(nondmg_aggregate, file.path(OUT_DIR, "tea_nondamaged_aggregate_summary.tsv"), sep = "\t")

fwrite(module_comp, file.path(OUT_DIR, "tea_nondamaged_vs_module_comparison.tsv"), sep = "\t")
fwrite(nondmg_vs_module_by_index, file.path(OUT_DIR, "tea_nondamaged_vs_module_by_index.tsv"), sep = "\t")
fwrite(brown_vs_nondmg, file.path(OUT_DIR, "tea_brown_vs_nondamaged_comparison.tsv"), sep = "\t")
fwrite(decision_notes, file.path(OUT_DIR, "04b_tea_variability_and_nondamaged_decision_notes.tsv"), sep = "\t")

write_run_metadata(
  file.path(OUT_DIR, "04b_tea_variability_and_nondamaged_run_metadata.tsv"),
  "04b_tea_variability_and_nondamaged.R",
  extra = list(
    max_age_kyr = PARAMS$main_max_age_kyr,
    n_module_variability_rows = nrow(module_variability),
    n_state_module_rows = nrow(state_module_variability),
    n_temporal_rows = nrow(temporal_summary),
    n_nondamaged_taxa = uniqueN(nondmg_taxon_annot$taxon),
    n_nondamaged_samples = nrow(nondmg_sample),
    nondamaged_light_filter = "none",
    comparison_focus_module = "brown",
    interpretation_scope = "descriptive_teapotential_not_activity"
  )
)
write_session_info(file.path(OUT_DIR, "04b_tea_variability_and_nondamaged_sessionInfo.txt"))

log_msg("04b_tea_variability_and_nondamaged complete.")
