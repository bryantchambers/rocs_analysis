library(lvplot)
library(showtext)
library(tidyverse)
library(DescTools)
library(ggh4x)
library(gghighlight)
library(data.table)
library(dtplyr)

# httpgd::hgd()
# httpgd::hgd_browse()

source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")
setwd("/projects/caeg/people/ngm902/apps/repos/rocs")


# Let's load the cdata
cdata <- read_tsv(file = "data/metadata_v5.tsv")



############################
## Species level data ######
############################

spdata <- read_tsv("taxonomic-profiling/eukaryotes/processed_data_june2026_incremental/merged_stats_metaDMG_species_level.tsv") %>%
  rename(label = library_id) %>%
  left_join(cdata) %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(total_n_reads_from_mD >= 100, mean_rlen > 35) %>%
  select(
    label, core, y_bp, depth_in_core_cm,
    # tax
    taxid, kingdom = kingdom_unicorn, phylum = phylum_unicorn , class = class_unicorn, order = order_unicorn, family = family_unicorn, genus = genus_unicorn, species = species.y_unicorn,
    # mapping stats
    n_references, num_alns, total_reads_from_unicorn, n_reads = total_n_reads_from_mD, mean_rlen, breath_ratio, duplicity,
    # damage
    Zfit = Zfit_new, A_b, q_b, c_b, phi_b, rho_c, rho_c_perm, rho_c_perm_pval)

# Statistics for all cores but ST5
spdata_stats <- spdata %>%
  filter(core != "ST5") %>%
  mutate(is_dmg = ifelse(Zfit >= 2 & A_b >= 0.1, "Damaged", "Non-damaged")) %>%
  group_by(kingdom, phylum, class, order, family, genus, species) %>%
  summarise(
    total_reads = sum(n_reads, na.rm = TRUE), 
    mean_rlen = mean(mean_rlen, na.rm = TRUE), 
    # mean_rlen = weighted.mean(mean_rlen, n_reads, na.rm = TRUE),
    n_labels = n_distinct(label), 
    mean_A_b = mean(A_b, na.rm = TRUE),
    # mean_A_b = weighted.mean(A_b, n_reads, na.rm = TRUE),
    median_A_b = median(A_b, na.rm = TRUE),
    sd_A_b = sd(A_b, na.rm = TRUE),
    mean_rho_c = mean(rho_c, na.rm = TRUE),
    # mean_rho_c = weighted.mean(rho_c, n_reads, na.rm = TRUE),
    mean_duplicity = mean(duplicity, na.rm = TRUE),
    # mean_duplicity = weighted.mean(duplicity, n_reads, na.rm = TRUE),
    mean_rho_c_perm = mean(rho_c_perm, na.rm = TRUE),
    # mean_rho_c_perm = weighted.mean(rho_c_perm, n_reads, na.rm = TRUE),
    mean_rho_c_perm_pval = mean(rho_c_perm_pval, na.rm = TRUE),
    # mean_rho_c_perm_pval = weighted.mean(rho_c_perm_pval, n_reads, na.rm = TRUE),
    proportion_damaged = sum(is_dmg == "Damaged") / n(),
    .groups = "drop")

filtered_species <- spdata_stats %>% 
  filter(
    n_labels > 5,
    total_reads > 1000,
    mean_A_b > 0.05,
    proportion_damaged > 0.5,
    mean_rho_c_perm > 0.5,
    mean_duplicity < 0.96
    )

## Statistics for ST5 only
spdata_stats_st5 <- spdata %>%
  filter(core == "ST5") %>%
  mutate(is_dmg = ifelse(Zfit >= 2 & A_b >= 0.1, "Damaged", "Non-damaged")) %>%
  group_by(kingdom, phylum, class, order, family, genus, species) %>%
  summarise(
    total_reads = sum(n_reads, na.rm = TRUE), 
    mean_rlen = mean(mean_rlen, na.rm = TRUE), 
    # mean_rlen = weighted.mean(mean_rlen, n_reads, na.rm = TRUE),
    n_labels = n_distinct(label), 
    mean_A_b = mean(A_b, na.rm = TRUE),
    # mean_A_b = weighted.mean(A_b, n_reads, na.rm = TRUE),
    median_A_b = median(A_b, na.rm = TRUE),
    sd_A_b = sd(A_b, na.rm = TRUE),
    mean_rho_c = mean(rho_c, na.rm = TRUE),
    # mean_rho_c = weighted.mean(rho_c, n_reads, na.rm = TRUE),
    mean_duplicity = mean(duplicity, na.rm = TRUE),
    # mean_duplicity = weighted.mean(duplicity, n_reads, na.rm = TRUE),
    mean_rho_c_perm = mean(rho_c_perm, na.rm = TRUE),
    # mean_rho_c_perm = weighted.mean(rho_c_perm, n_reads, na.rm = TRUE),
    mean_rho_c_perm_pval = mean(rho_c_perm_pval, na.rm = TRUE),
    # mean_rho_c_perm_pval = weighted.mean(rho_c_perm_pval, n_reads, na.rm = TRUE),
    proportion_damaged = sum(is_dmg == "Damaged") / n(),
    .groups = "drop")

filtered_species_st5 <- spdata_stats_st5 %>% 
  filter(
    n_labels > 5,
    total_reads > 1000,
    mean_A_b > 0.05,
    proportion_damaged > 0.5,
    mean_rho_c_perm > 0.5,
    mean_duplicity < 0.96
    )


## Apply filterings
spdata <- spdata %>%
    mutate(
      is_dmg_local = ifelse(Zfit > 2 & A_b >= 0.1, "Damaged", "Non-damaged"),
      # is_dmg = ifelse(species %in% filtered_species$species, "Damaged", "Non-damaged")
      is_dmg = case_when(
        core != "ST5" & species %in% filtered_species$species ~ "Damaged",
        core == "ST5" & species %in% filtered_species_st5$species ~ "Damaged",
        TRUE ~ "Non-damaged"
        )
      ) 

spdata %>%
  write_tsv(., "/projects/caeg/people/ngm902/apps/repos/rocs/results/eukaryotes/taxonomy/species_level_data.tsv")

spdata %>%
  filter(is_dmg == "Damaged") %>%
  write_tsv(., "/projects/caeg/people/ngm902/apps/repos/rocs/results/eukaryotes/taxonomy/species_level_data_dmg.tsv")


############################
## Genus level data ########
############################

gendata <- read_tsv("taxonomic-profiling/eukaryotes/processed_data_june2026_incremental/merged_stats_metaDMG_genus_level.tsv") %>%
  rename(label = library_id) %>%
  left_join(cdata) %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(total_n_reads_from_mD >= 100 & mean_rlen > 35) %>%
  select(
    label, core, y_bp,
    # tax
    taxid, kingdom = kingdom_unicorn, phylum = phylum_unicorn , class = class_unicorn, order = order_unicorn, family = family_unicorn, genus = genus.x_unicorn,
    # mapping stats
    n_references, num_alns, total_reads_from_unicorn, n_reads = total_n_reads_from_mD, mean_rlen, breath_ratio, duplicity,
    # damage
    Zfit = Zfit_new, A_b, q_b, c_b, phi_b, rho_c, rho_c_perm, rho_c_perm_pval)

# Statistics for all cores but ST5
gendata_stats <- gendata %>%
  filter(core != "ST5") %>%
  mutate(is_dmg = ifelse(Zfit >= 2 & A_b >= 0.1, "Damaged", "Non-damaged")) %>%
  group_by(kingdom, phylum, class, order, family, genus) %>%
  summarise(
    total_reads = sum(n_reads, na.rm = TRUE), 
    mean_rlen = mean(mean_rlen, na.rm = TRUE), 
    # mean_rlen = weighted.mean(mean_rlen, n_reads, na.rm = TRUE),
    n_labels = n_distinct(label), 
    mean_A_b = mean(A_b, na.rm = TRUE),
    # mean_A_b = weighted.mean(A_b, n_reads, na.rm = TRUE),
    median_A_b = median(A_b, na.rm = TRUE),
    sd_A_b = sd(A_b, na.rm = TRUE),
    mean_rho_c = mean(rho_c, na.rm = TRUE),
    # mean_rho_c = weighted.mean(rho_c, n_reads, na.rm = TRUE),
    mean_duplicity = mean(duplicity, na.rm = TRUE),
    # mean_duplicity = weighted.mean(duplicity, n_reads, na.rm = TRUE),
    mean_rho_c_perm = mean(rho_c_perm, na.rm = TRUE),
    # mean_rho_c_perm = weighted.mean(rho_c_perm, n_reads, na.rm = TRUE),
    mean_rho_c_perm_pval = mean(rho_c_perm_pval, na.rm = TRUE),
    # mean_rho_c_perm_pval = weighted.mean(rho_c_perm_pval, n_reads, na.rm = TRUE),
    proportion_damaged = sum(is_dmg == "Damaged") / n(),
    .groups = "drop")

filtered_genera <- gendata_stats %>%
  filter(
    n_labels > 5,
    total_reads > 1000,
    mean_A_b > 0.05,
    proportion_damaged > 0.5,
    mean_rho_c_perm > 0.5,
    mean_duplicity < 0.96
    )

# Statistics for ST5
gendata_stats_st5 <- gendata %>%
  filter(core == "ST5") %>%
  mutate(is_dmg = ifelse(Zfit >= 2 & A_b >= 0.1, "Damaged", "Non-damaged")) %>%
  group_by(kingdom, phylum, class, order, family, genus) %>%
  summarise(
    total_reads = sum(n_reads, na.rm = TRUE), 
    mean_rlen = mean(mean_rlen, na.rm = TRUE), 
    # mean_rlen = weighted.mean(mean_rlen, n_reads, na.rm = TRUE),
    n_labels = n_distinct(label), 
    mean_A_b = mean(A_b, na.rm = TRUE),
    # mean_A_b = weighted.mean(A_b, n_reads, na.rm = TRUE),
    median_A_b = median(A_b, na.rm = TRUE),
    sd_A_b = sd(A_b, na.rm = TRUE),
    mean_rho_c = mean(rho_c, na.rm = TRUE),
    # mean_rho_c = weighted.mean(rho_c, n_reads, na.rm = TRUE),
    mean_duplicity = mean(duplicity, na.rm = TRUE),
    # mean_duplicity = weighted.mean(duplicity, n_reads, na.rm = TRUE),
    mean_rho_c_perm = mean(rho_c_perm, na.rm = TRUE),
    # mean_rho_c_perm = weighted.mean(rho_c_perm, n_reads, na.rm = TRUE),
    mean_rho_c_perm_pval = mean(rho_c_perm_pval, na.rm = TRUE),
    # mean_rho_c_perm_pval = weighted.mean(rho_c_perm_pval, n_reads, na.rm = TRUE),
    proportion_damaged = sum(is_dmg == "Damaged") / n(),
    .groups = "drop")

filtered_genera_st5 <- gendata_stats_st5 %>%
  filter(
    n_labels > 5,
    total_reads > 1000,
    mean_A_b > 0.05,
    proportion_damaged > 0.5,
    mean_rho_c_perm > 0.5,
    mean_duplicity < 0.96
    )

## Apply filterings
gendata <- gendata %>%
    mutate(
      is_dmg_local = ifelse(Zfit > 2 & A_b >= 0.1, "Damaged", "Non-damaged"),
      is_dmg = case_when(
          core != "ST5" & genus %in% filtered_genera$genus ~ "Damaged",
          core == "ST5" & genus %in% filtered_genera_st5$genus ~ "Damaged",
          TRUE ~ "Non-damaged"
          )
        )  

gendata %>%
  write_tsv(., "/projects/caeg/people/ngm902/apps/repos/rocs/results/eukaryotes/taxonomy/genus_level_data.tsv")

gendata %>%
  filter(is_dmg == "Damaged") %>%
  write_tsv(., "/projects/caeg/people/ngm902/apps/repos/rocs/results/eukaryotes/taxonomy/genus_level_data_dmg.tsv")


############################
## Family level data ########
############################

famdata <- read_tsv("taxonomic-profiling/eukaryotes/processed_data_june2026_incremental/merged_stats_metaDMG_family_level.tsv") %>%
  rename(label = library_id) %>%
  left_join(cdata) %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(total_n_reads_from_mD >= 100 & mean_rlen > 35) %>%
  select(
    label, core, y_bp,
    # tax
    taxid, kingdom = kingdom_unicorn, phylum = phylum_unicorn , class = class_unicorn, order = order_unicorn, family = family.x_unicorn,
    # mapping stats
    n_references, num_alns, total_reads_from_unicorn, n_reads = total_n_reads_from_mD, mean_rlen, breath_ratio, duplicity,
    # damage
    Zfit = Zfit_new, A_b, q_b, c_b, phi_b, rho_c, rho_c_perm, rho_c_perm_pval)

# Statistics for all cores but ST5
famdata_stats <- famdata %>%
  filter(core != "ST5") %>%
  mutate(is_dmg = ifelse(Zfit >= 2 & A_b >= 0.1, "Damaged", "Non-damaged")) %>%
  group_by(kingdom, phylum, class, order, family) %>%
  summarise(
    total_reads = sum(n_reads, na.rm = TRUE), 
    mean_rlen = mean(mean_rlen, na.rm = TRUE), 
    # mean_rlen = weighted.mean(mean_rlen, n_reads, na.rm = TRUE),
    n_labels = n_distinct(label), 
    mean_A_b = mean(A_b, na.rm = TRUE),
    # mean_A_b = weighted.mean(A_b, n_reads, na.rm = TRUE),
    median_A_b = median(A_b, na.rm = TRUE),
    sd_A_b = sd(A_b, na.rm = TRUE),
    mean_rho_c = mean(rho_c, na.rm = TRUE),
    # mean_rho_c = weighted.mean(rho_c, n_reads, na.rm = TRUE),
    mean_duplicity = mean(duplicity, na.rm = TRUE),
    # mean_duplicity = weighted.mean(duplicity, n_reads, na.rm = TRUE),
    mean_rho_c_perm = mean(rho_c_perm, na.rm = TRUE),
    # mean_rho_c_perm = weighted.mean(rho_c_perm, n_reads, na.rm = TRUE),
    mean_rho_c_perm_pval = mean(rho_c_perm_pval, na.rm = TRUE),
    # mean_rho_c_perm_pval = weighted.mean(rho_c_perm_pval, n_reads, na.rm = TRUE),
    proportion_damaged = sum(is_dmg == "Damaged") / n(),
    .groups = "drop")

filtered_family <- famdata_stats %>%
  filter(
    n_labels > 5,
    total_reads > 1000,
    mean_A_b > 0.05,
    proportion_damaged > 0.5,
    mean_rho_c_perm > 0.5,
    mean_duplicity < 0.96
    )

# Statistics for ST5
famdata_stats_st5 <- famdata %>%
  filter(core == "ST5") %>%
  mutate(is_dmg = ifelse(Zfit >= 2 & A_b >= 0.1, "Damaged", "Non-damaged")) %>%
  group_by(kingdom, phylum, class, order, family) %>%
  summarise(
    total_reads = sum(n_reads, na.rm = TRUE), 
    mean_rlen = mean(mean_rlen, na.rm = TRUE), 
    # mean_rlen = weighted.mean(mean_rlen, n_reads, na.rm = TRUE),
    n_labels = n_distinct(label), 
    mean_A_b = mean(A_b, na.rm = TRUE),
    # mean_A_b = weighted.mean(A_b, n_reads, na.rm = TRUE),
    median_A_b = median(A_b, na.rm = TRUE),
    sd_A_b = sd(A_b, na.rm = TRUE),
    mean_rho_c = mean(rho_c, na.rm = TRUE),
    # mean_rho_c = weighted.mean(rho_c, n_reads, na.rm = TRUE),
    mean_duplicity = mean(duplicity, na.rm = TRUE),
    # mean_duplicity = weighted.mean(duplicity, n_reads, na.rm = TRUE),
    mean_rho_c_perm = mean(rho_c_perm, na.rm = TRUE),
    # mean_rho_c_perm = weighted.mean(rho_c_perm, n_reads, na.rm = TRUE),
    mean_rho_c_perm_pval = mean(rho_c_perm_pval, na.rm = TRUE),
    # mean_rho_c_perm_pval = weighted.mean(rho_c_perm_pval, n_reads, na.rm = TRUE),
    proportion_damaged = sum(is_dmg == "Damaged") / n(),
    .groups = "drop")

filtered_family_st5 <- famdata_stats_st5 %>%
  filter(
    n_labels > 5,
    total_reads > 1000,
    mean_A_b > 0.05,
    proportion_damaged > 0.5,
    mean_rho_c_perm > 0.5,
    mean_duplicity < 0.96
    )

famdata <- famdata %>%
    mutate(
      is_dmg_local = ifelse(Zfit > 2 & A_b >= 0.1, "Damaged", "Non-damaged"),
      is_dmg = case_when(
          core != "ST5" & family %in% filtered_family$family ~ "Damaged",
          core == "ST5" & family %in% filtered_family_st5$family ~ "Damaged",
          TRUE ~ "Non-damaged"
          )
        )  

famdata %>%
  write_tsv(., "/projects/caeg/people/ngm902/apps/repos/rocs/results/eukaryotes/taxonomy/family_level_data.tsv")

famdata %>%
  filter(is_dmg == "Damaged") %>%
  write_tsv(., "/projects/caeg/people/ngm902/apps/repos/rocs/results/eukaryotes/taxonomy/family_level_data_dmg.tsv")
