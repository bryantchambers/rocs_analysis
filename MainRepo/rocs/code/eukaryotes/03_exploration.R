#!/usr/bin/env Rscript

# =============================================================================
# Eukaryote exploratory analysis (publication-oriented, conservative framing)
#
# Notes:
# - Damage-supported detections are treated as candidate ancient signals.
# - Detection is not equivalent to local abundance (especially for vertebrates).
# - Compositional and uneven-coverage caveats apply to all read-count summaries.
# =============================================================================

setwd("/maps/projects/caeg/people/ngm902/apps/repos/rocs")
repo_root <- getwd()

suppressPackageStartupMessages({
  library(tidyverse)
})

analysis_notes <- c(
  "Damage-supported detections are candidate ancient signals, not proof of local abundance.",
  "Read counts are compositional; relative abundance summaries are exploratory only.",
  "Coverage and richness comparisons can be biased by sequencing depth and preservation."
)

required_files <- c(
  "data/metadata_v5.tsv",
  "results/eukaryotes/taxonomy/species_level_data.tsv",
  "results/eukaryotes/taxonomy/genus_level_data.tsv",
  "results/eukaryotes/taxonomy/family_level_data.tsv"
)

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop(
    "Missing required input files:\n - ",
    paste(missing_required, collapse = "\n - "),
    call. = FALSE
  )
}

# ------------------------------
# Core analysis settings
# ------------------------------
core_levels <- c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")
taxonomy_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# ------------------------------
# Helper functions
# ------------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

alpha_metrics <- function(counts) {
  counts <- safe_numeric(counts)
  counts <- counts[is.finite(counts) & counts > 0]
  n <- sum(counts)
  if (!is.finite(n) || n <= 0) {
    return(tibble(
      observed_richness = 0,
      Shannon = NA_real_,
      InvSimpson = NA_real_
    ))
  }

  p <- counts / n
  shannon <- -sum(p * log(p))
  inv_simpson <- ifelse(sum(p^2) > 0, 1 / sum(p^2), NA_real_)

  tibble(
    observed_richness = length(counts),
    Shannon = shannon,
    InvSimpson = inv_simpson
  )
}

collapse_taxonomy_text <- function(...) {
  vals <- c(...)
  vals <- as.character(vals)
  vals[is.na(vals)] <- ""
  vals <- stringr::str_squish(tolower(vals))
  vals <- vals[vals != ""]
  if (length(vals) == 0) return("")
  paste(vals, collapse = " ; ")
}

classify_ecological_group <- function(df) {
  tx <- df %>% select(any_of(taxonomy_cols))
  tax_text <- pmap_chr(
    tx,
    collapse_taxonomy_text
  )

  phytoplankton_pattern <- paste0(
    "bacillariophyta|dinophyta|haptophyta|chlorophyta|cryptophyta|ochrophyta|prasinophy|",
    "phaeocyst|emiliania|coccolith|chaetoceros|thalassiosira|skeletonema|triparmaceae|noelaerhabdaceae|",
    "aureococcus|pelagomonas|pelagomonadaceae|teleaulax|florenciella|octactis|dictyochaceae|geminigeraceae|symbiodiniaceae"
  )
  macroalgae_pattern <- paste0(
    "phaeophyceae|laminariales|laminariaceae|saccharina|kelp"
  )
  zooplankton_pattern <- paste0(
    "copepoda|hexanauplia|calanoida|euphausiacea|krill|appendicularia|chaetognatha|",
    "foraminifera|radiolaria|ostracoda|siphonophora|pteropod|salp|cladocera|",
    "ctenophora|tentaculata|cydippida|lobata|mertensia|mertensiidae|bolinopsis|bolinopsidae|",
    "collozoum|collodaria|polycystinea|sphaerozoidae|nausithoe|nausithoidae"
  )
  fish_pattern <- paste0(
    "actinopteri|chondrichthyes|myctophiformes|gadiformes|clupeiformes|perciformes|salmoniformes|",
    "myctophidae|gadidae|clupeidae|scombridae"
  )
  marine_mammal_pattern <- paste0(
    "cetacea|mysticeti|odontoceti|phocidae|otariidae|odobenidae|delphinidae|balaenidae|physeteridae|",
    "balaenoptera|megaptera|eschrichtius|hyperoodon|mesoplodon|balaenopteridae|ziphiidae|monodontidae|phocoenidae|eschrichtiidae"
  )
  benthic_pattern <- paste0(
    "bivalvia|gastropoda|cephalopoda|polychaeta|annelida|echinodermata|asteroidea|ophiuroidea|",
    "holothuroidea|echinoidea|porifera|bryozoa|ascidiacea|actiniaria|malacostraca|decapoda|amphipoda|isopoda|",
    "chirona|semibalanus|balanus|balanidae|chrysogorgia|chrysogorgiidae|coralliidae"
  )
  non_marine_ambig_pattern <- paste0(
    "primates|hominidae|aves|insecta|arachnida|embryophyta|streptophyta|fungi|amphibia|reptilia|",
    "\\bmus\\b|\\bgallus\\b|\\bcanis\\b|\\bfelis\\b|\\bsus\\b|mustelidae|paramecium|parameciidae|uncultured eukaryote"
  )

  case_when(
    str_detect(tax_text, phytoplankton_pattern) ~ "Phytoplankton",

    str_detect(tax_text, macroalgae_pattern) ~ "Macroalgae",

    str_detect(tax_text, zooplankton_pattern) ~ "Zooplankton",

    str_detect(tax_text, fish_pattern) ~ "Fish",

    str_detect(tax_text, marine_mammal_pattern) ~ "Marine mammals",

    str_detect(tax_text, benthic_pattern) ~ "Benthic invertebrates",

    str_detect(tax_text, non_marine_ambig_pattern) ~ "Ambiguous/non-marine",

    TRUE ~ "Other"
  )
}

top_taxa_summary <- function(df, top_n = 20) {
  if (nrow(df) == 0) return(tibble())

  sample_totals <- df %>%
    group_by(label) %>%
    summarise(sample_total = sum(n_reads, na.rm = TRUE), .groups = "drop")

  tax_sum <- df %>%
    left_join(sample_totals, by = "label") %>%
    mutate(
      rel_abundance = if_else(sample_total > 0, n_reads / sample_total, NA_real_),
      present = n_reads > 0
    ) %>%
    group_by(taxon, ecological_group) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      n_samples = sum(present, na.rm = TRUE),
      prevalence = n_samples / n_distinct(label),
      median_rel_abundance_when_present = safe_median(rel_abundance[present]),
      max_rel_abundance = safe_max(rel_abundance),
      .groups = "drop"
    ) %>%
    arrange(desc(total_reads), desc(prevalence))

  tax_sum %>%
    mutate(rank_order = row_number()) %>%
    filter(rank_order <= top_n) %>%
    select(-rank_order)
}

load_and_standardize_rank <- function(path, rank, taxon_col, metadata_tbl) {
  raw <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)

  if (!"label" %in% names(raw)) stop("Missing `label` in: ", path, call. = FALSE)
  if (!"n_reads" %in% names(raw)) stop("Missing `n_reads` in: ", path, call. = FALSE)
  if (!taxon_col %in% names(raw)) stop("Missing rank column `", taxon_col, "` in: ", path, call. = FALSE)

  missing_tax_cols <- setdiff(taxonomy_cols, names(raw))
  for (col in missing_tax_cols) raw[[col]] <- NA_character_

  if (!"is_dmg" %in% names(raw)) raw$is_dmg <- NA_character_
  if (!"core" %in% names(raw)) raw$core <- NA_character_
  if (!"y_bp" %in% names(raw)) raw$y_bp <- NA_real_
  if (!"depth_in_core_cm" %in% names(raw)) raw$depth_in_core_cm <- NA_real_

  out <- raw %>%
    mutate(
      n_reads = safe_numeric(n_reads),
      y_bp = safe_numeric(y_bp),
      depth_in_core_cm = safe_numeric(depth_in_core_cm),
      rank = rank,
      taxon = as.character(.data[[taxon_col]]),
      taxon = if_else(is.na(taxon) | taxon == "", "Unclassified", taxon),
      is_dmg = as.character(is_dmg)
    ) %>%
    left_join(
      metadata_tbl %>%
        select(label, core_meta = core, y_bp_meta = y_bp, depth_meta = depth_in_core_cm),
      by = "label"
    ) %>%
    mutate(
      core = coalesce(as.character(core), as.character(core_meta)),
      y_bp = coalesce(y_bp, safe_numeric(y_bp_meta)),
      depth_in_core_cm = coalesce(depth_in_core_cm, safe_numeric(depth_meta))
    ) %>%
    select(-core_meta, -y_bp_meta, -depth_meta) %>%
    filter(core %in% core_levels) %>%
    mutate(
      core = factor(core, levels = core_levels)
    )

  out$ecological_group <- classify_ecological_group(out)

  list(
    full = out,
    dmg = out %>% filter(is_dmg == "Damaged")
  )
}

build_rank_summaries <- function(df, rank_name) {
  if (nrow(df) == 0) {
    return(list(
      sample_summary = tibble(),
      taxon_summary = tibble(),
      group_summary_dmg = tibble(),
      sample_vs_depth_summary = tibble()
    ))
  }

  sample_summary <- df %>%
    group_by(label, core, y_bp, depth_in_core_cm) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      n_taxa = n_distinct(taxon[n_reads > 0]),
      count_vector = list(n_reads[n_reads > 0]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(alpha = list(alpha_metrics(unlist(count_vector)))) %>%
    tidyr::unnest_wider(alpha) %>%
    ungroup() %>%
    select(-count_vector) %>%
    mutate(rank = rank_name)

  sample_totals <- df %>%
    group_by(label) %>%
    summarise(sample_total = sum(n_reads, na.rm = TRUE), .groups = "drop")

  taxon_summary <- df %>%
    left_join(sample_totals, by = "label") %>%
    mutate(
      present = n_reads > 0,
      rel_abundance = if_else(sample_total > 0, n_reads / sample_total, NA_real_)
    ) %>%
    group_by(taxon, ecological_group, across(any_of(taxonomy_cols))) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      n_samples = sum(present, na.rm = TRUE),
      prevalence = n_samples / n_distinct(label),
      median_rel_abundance_when_present = safe_median(rel_abundance[present]),
      max_rel_abundance = safe_max(rel_abundance),
      .groups = "drop"
    ) %>%
    arrange(desc(total_reads), desc(prevalence)) %>%
    mutate(rank = rank_name)

  group_summary <- df %>%
    group_by(label) %>%
    mutate(sample_total = sum(n_reads, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(label, core, y_bp, depth_in_core_cm, ecological_group) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      n_taxa = n_distinct(taxon[n_reads > 0]),
      sample_total = first(sample_total),
      .groups = "drop"
    ) %>%
    mutate(relative_abundance = if_else(sample_total > 0, total_reads / sample_total, NA_real_)) %>%
    select(-sample_total) %>%
    mutate(rank = rank_name)

  sample_vs_depth_summary <- sample_summary %>%
    mutate(
      log10_total_reads = log10(total_reads + 1),
      log10_observed_richness = log10(observed_richness + 1)
    ) %>%
    select(rank, label, core, y_bp, depth_in_core_cm, total_reads, observed_richness,
           Shannon, InvSimpson, log10_total_reads, log10_observed_richness)

  list(
    sample_summary = sample_summary,
    taxon_summary = taxon_summary,
    group_summary_dmg = group_summary,
    sample_vs_depth_summary = sample_vs_depth_summary
  )
}

summarize_damaged_vs_all <- function(rank_data_list) {
  bind_rows(lapply(names(rank_data_list), function(rank_name) {
    full <- rank_data_list[[rank_name]]$full
    dmg <- rank_data_list[[rank_name]]$dmg

    full_reads <- sum(full$n_reads, na.rm = TRUE)
    dmg_reads <- sum(dmg$n_reads, na.rm = TRUE)
    full_taxa <- n_distinct(full$taxon[full$n_reads > 0])
    dmg_taxa <- n_distinct(dmg$taxon[dmg$n_reads > 0])

    tibble(
      rank = rank_name,
      n_labels_full = n_distinct(full$label),
      n_labels_dmg = n_distinct(dmg$label),
      total_reads_full = full_reads,
      total_reads_dmg = dmg_reads,
      frac_reads_dmg = ifelse(full_reads > 0, dmg_reads / full_reads, NA_real_),
      n_taxa_full = full_taxa,
      n_taxa_dmg = dmg_taxa,
      frac_taxa_dmg = ifelse(full_taxa > 0, dmg_taxa / full_taxa, NA_real_)
    )
  }))
}

build_ambiguous_sample_flags <- function(df, rel_abundance_flag = 0.2) {
  if (nrow(df) == 0) return(tibble())

  by_sample <- df %>%
    group_by(label, core, y_bp, depth_in_core_cm) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      ambiguous_reads = sum(n_reads[ecological_group == "Ambiguous/non-marine"], na.rm = TRUE),
      ambiguous_taxa = n_distinct(taxon[ecological_group == "Ambiguous/non-marine" & n_reads > 0]),
      .groups = "drop"
    ) %>%
    mutate(
      ambiguous_rel_abundance = if_else(total_reads > 0, ambiguous_reads / total_reads, NA_real_),
      high_ambiguous_burden = ambiguous_rel_abundance >= rel_abundance_flag
    )

  top_ambiguous <- df %>%
    filter(ecological_group == "Ambiguous/non-marine") %>%
    group_by(label, taxon) %>%
    summarise(taxon_reads = sum(n_reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(label) %>%
    arrange(desc(taxon_reads), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    rename(top_ambiguous_taxon = taxon, top_ambiguous_taxon_reads = taxon_reads)

  by_sample %>%
    left_join(top_ambiguous, by = "label") %>%
    arrange(desc(ambiguous_rel_abundance), desc(ambiguous_reads))
}

build_cross_rank_concordance <- function(sample_summary_dmg_df) {
  rank_wide <- sample_summary_dmg_df %>%
    select(label, core, y_bp, rank, observed_richness, total_reads, Shannon) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = rank,
      values_from = c(observed_richness, total_reads, Shannon),
      names_sep = "__"
    )

  pairwise <- tribble(
    ~rank1, ~rank2,
    "species", "genus",
    "species", "family",
    "genus", "family"
  )

  bind_rows(lapply(seq_len(nrow(pairwise)), function(i) {
    r1 <- pairwise$rank1[i]
    r2 <- pairwise$rank2[i]

    tibble(
      rank1 = r1,
      rank2 = r2,
      rho_observed_richness = suppressWarnings(cor(
        rank_wide[[paste0("observed_richness__", r1)]],
        rank_wide[[paste0("observed_richness__", r2)]],
        method = "spearman",
        use = "pairwise.complete.obs"
      )),
      rho_total_reads = suppressWarnings(cor(
        rank_wide[[paste0("total_reads__", r1)]],
        rank_wide[[paste0("total_reads__", r2)]],
        method = "spearman",
        use = "pairwise.complete.obs"
      )),
      rho_shannon = suppressWarnings(cor(
        rank_wide[[paste0("Shannon__", r1)]],
        rank_wide[[paste0("Shannon__", r2)]],
        method = "spearman",
        use = "pairwise.complete.obs"
      ))
    )
  }))
}

build_suspicious_taxa_summary <- function(df, top_n = 25) {
  if (nrow(df) == 0) return(tibble())

  sample_total <- n_distinct(df$label)
  df %>%
    filter(ecological_group == "Ambiguous/non-marine") %>%
    group_by(taxon, across(any_of(taxonomy_cols))) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      n_samples = n_distinct(label[n_reads > 0]),
      prevalence = ifelse(sample_total > 0, n_samples / sample_total, NA_real_),
      max_reads_single_sample = safe_max(n_reads),
      .groups = "drop"
    ) %>%
    arrange(desc(total_reads), desc(n_samples)) %>%
    slice_head(n = top_n)
}

# ------------------------------
# Data loading
# ------------------------------
metadata_raw <- readr::read_tsv("data/metadata_v5.tsv", show_col_types = FALSE, progress = FALSE)
metadata_used <- metadata_raw %>%
  mutate(
    y_bp = safe_numeric(y_bp),
    depth_in_core_cm = safe_numeric(depth_in_core_cm),
    core = as.character(core)
  ) %>%
  filter(core %in% core_levels) %>%
  mutate(core = factor(core, levels = core_levels))

species_data <- load_and_standardize_rank(
  path = "results/eukaryotes/taxonomy/species_level_data.tsv",
  rank = "species",
  taxon_col = "species",
  metadata_tbl = metadata_used
)

genus_data <- load_and_standardize_rank(
  path = "results/eukaryotes/taxonomy/genus_level_data.tsv",
  rank = "genus",
  taxon_col = "genus",
  metadata_tbl = metadata_used
)

family_data <- load_and_standardize_rank(
  path = "results/eukaryotes/taxonomy/family_level_data.tsv",
  rank = "family",
  taxon_col = "family",
  metadata_tbl = metadata_used
)

rank_data <- list(species = species_data, genus = genus_data, family = family_data)

# ------------------------------
# Rank-level summaries
# ------------------------------
species_summaries_full <- build_rank_summaries(rank_data$species$full, "species")
species_summaries_dmg <- build_rank_summaries(rank_data$species$dmg, "species")

species_summary_full <- species_summaries_full$sample_summary
species_summary_dmg <- species_summaries_dmg$sample_summary
species_taxon_full <- species_summaries_full$taxon_summary
species_taxon_dmg <- species_summaries_dmg$taxon_summary
species_group_dmg <- species_summaries_dmg$group_summary_dmg
species_depth_qc <- species_summaries_full$sample_vs_depth_summary

genus_summaries_full <- build_rank_summaries(rank_data$genus$full, "genus")
genus_summaries_dmg <- build_rank_summaries(rank_data$genus$dmg, "genus")

family_summaries_full <- build_rank_summaries(rank_data$family$full, "family")
family_summaries_dmg <- build_rank_summaries(rank_data$family$dmg, "family")

sample_summary_full <- bind_rows(
  species_summary_full,
  genus_summaries_full$sample_summary,
  family_summaries_full$sample_summary
)

sample_summary_dmg <- bind_rows(
  species_summary_dmg,
  genus_summaries_dmg$sample_summary,
  family_summaries_dmg$sample_summary
)

taxon_summary_full <- bind_rows(
  species_taxon_full,
  genus_summaries_full$taxon_summary,
  family_summaries_full$taxon_summary
)

taxon_summary_dmg <- bind_rows(
  species_taxon_dmg,
  genus_summaries_dmg$taxon_summary,
  family_summaries_dmg$taxon_summary
)

group_summary_dmg <- bind_rows(
  species_group_dmg,
  genus_summaries_dmg$group_summary_dmg,
  family_summaries_dmg$group_summary_dmg
)

sample_vs_depth_summary <- bind_rows(
  species_depth_qc,
  genus_summaries_full$sample_vs_depth_summary,
  family_summaries_full$sample_vs_depth_summary
)

damaged_vs_all_summary <- summarize_damaged_vs_all(rank_data)

cross_rank_concordance <- build_cross_rank_concordance(sample_summary_dmg)

suspicious_taxa_top_genus_dmg <- build_suspicious_taxa_summary(rank_data$genus$dmg, top_n = 250)
suspicious_taxa_top_species_dmg <- build_suspicious_taxa_summary(rank_data$species$dmg, top_n = 250)

ambiguous_sample_flags <- build_ambiguous_sample_flags(rank_data$genus$dmg, rel_abundance_flag = 0.2)

# Extra QC summary table from full data
read_richness_relationship <- sample_vs_depth_summary %>%
  group_by(rank, core) %>%
  summarise(
    n_samples = n(),
    cor_reads_richness = suppressWarnings(cor(total_reads, observed_richness, method = "spearman", use = "pairwise.complete.obs")),
    cor_reads_shannon = suppressWarnings(cor(total_reads, Shannon, method = "spearman", use = "pairwise.complete.obs")),
    .groups = "drop"
  )

depth_confounding_flags <- read_richness_relationship %>%
  mutate(
    abs_cor_reads_richness = abs(cor_reads_richness),
    abs_cor_reads_shannon = abs(cor_reads_shannon),
    strong_depth_richness_coupling = abs_cor_reads_richness >= 0.7,
    strong_depth_shannon_coupling = abs_cor_reads_shannon >= 0.7
  ) %>%
  arrange(desc(abs_cor_reads_richness))

if (nrow(depth_confounding_flags) > 0 && any(depth_confounding_flags$strong_depth_richness_coupling, na.rm = TRUE)) {
  flagged <- depth_confounding_flags %>%
    filter(strong_depth_richness_coupling) %>%
    transmute(id = paste0(rank, "@", as.character(core), " (rho=", signif(cor_reads_richness, 3), ")")) %>%
    pull(id)

  analysis_notes <- c(
    analysis_notes,
    paste0(
      "Strong sequencing-depth coupling with observed richness in: ",
      paste(flagged, collapse = "; "),
      ". Richness trends should be interpreted as coverage-sensitive."
    )
  )
}

if (nrow(ambiguous_sample_flags) > 0) {
  n_high_amb <- sum(ambiguous_sample_flags$high_ambiguous_burden, na.rm = TRUE)
  analysis_notes <- c(
    analysis_notes,
    paste0(
      "Samples with high ambiguous/non-marine burden (>=20% genus damaged reads): ",
      n_high_amb, " of ", nrow(ambiguous_sample_flags),
      "."
    )
  )
}

# ------------------------------
# Organize outputs
# ------------------------------
species_exploration <- list(
  data = rank_data$species,
  summaries = list(
    sample_summary_full = species_summary_full,
    sample_summary_dmg = species_summary_dmg,
    taxon_summary_full = species_taxon_full,
    taxon_summary_dmg = species_taxon_dmg,
    group_summary_dmg = species_group_dmg,
    sample_vs_depth_summary = species_depth_qc
  ),
  top_taxa = list(
    full = top_taxa_summary( rank_data$species$full, top_n = 25),
    dmg = top_taxa_summary(rank_data$species$dmg, top_n = 25)
  )
)
saveRDS(species_exploration, file = "results/eukaryotes/taxonomy/species_exploration.rds")

genus_exploration <- list(
  data = rank_data$genus,
  summaries = list(
    sample_summary_full = genus_summaries_full$sample_summary,
    sample_summary_dmg = genus_summaries_dmg$sample_summary,
    taxon_summary_full = genus_summaries_full$taxon_summary,
    taxon_summary_dmg = genus_summaries_dmg$taxon_summary,
    group_summary_dmg = genus_summaries_dmg$group_summary_dmg,
    sample_vs_depth_summary = genus_summaries_full$sample_vs_depth_summary
  ),
  top_taxa = list(
    full = top_taxa_summary(rank_data$genus$full, top_n = 25),
    dmg = top_taxa_summary(rank_data$genus$dmg, top_n = 25)
  )
)
saveRDS(genus_exploration, file = "results/eukaryotes/taxonomy/genus_exploration.rds")

family_exploration <- list(
  data = rank_data$family,
  summaries = list(
    sample_summary_full = family_summaries_full$sample_summary,
    sample_summary_dmg = family_summaries_dmg$sample_summary,
    taxon_summary_full = family_summaries_full$taxon_summary,
    taxon_summary_dmg = family_summaries_dmg$taxon_summary,
    group_summary_dmg = family_summaries_dmg$group_summary_dmg,
    sample_vs_depth_summary = family_summaries_full$sample_vs_depth_summary
  ),
  top_taxa = list(
    full = top_taxa_summary(rank_data$family$full, top_n = 25),
    dmg = top_taxa_summary(rank_data$family$dmg, top_n = 25)
  )
)
saveRDS(family_exploration, file = "results/eukaryotes/taxonomy/family_exploration.rds")

exploration_results <- list(
  metadata_used = metadata_used,
  settings = list(
    repo_root = repo_root,
    cores_used = core_levels,
    required_files = required_files
  ),
  data = list(
    species = rank_data$species,
    genus = rank_data$genus,
    family = rank_data$family
  ),
  summaries = list(
    sample_summary_full = sample_summary_full,
    sample_summary_dmg = sample_summary_dmg,
    taxon_summary_full = taxon_summary_full,
    taxon_summary_dmg = taxon_summary_dmg,
    group_summary_dmg = group_summary_dmg,
    sample_vs_depth_summary = sample_vs_depth_summary,
    read_richness_relationship = read_richness_relationship,
    depth_confounding_flags = depth_confounding_flags,
    damaged_vs_all_summary = damaged_vs_all_summary,
    cross_rank_concordance = cross_rank_concordance,
    suspicious_taxa_top_genus_dmg = suspicious_taxa_top_genus_dmg,
    ambiguous_sample_flags = ambiguous_sample_flags
  ),
  notes = unique(analysis_notes)
)

message("Exploration objects created in memory:")
message(" - exploration_results")
message(" - species_exploration")
message(" - genus_exploration")
message(" - family_exploration")




