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

spdata <- read_tsv("results/eukaryotes/taxonomy/species_level_data.tsv")

plotA <- spdata %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg_local = fct_relevel(is_dmg_local, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = A_b, fill = is_dmg_local, size = n_reads))+
      geom_point(color = "black", shape = 21, alpha = 0.5, stroke = 0.2)+
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
      facet_nested(.~core*is_dmg)+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))+
      labs(x = "Age (y bp)", y = "Damage (A_b)", color = "Local-based classification (A_b > 0.1 & Zfit >= 2)", size = "Number of reads")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/eukaryotes/preliminary/damage_depositional_classification_species.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()


plotA <- spdata %>%
    group_by(label, core, y_bp, is_dmg) %>%
    summarise(n_reads = sum(n_reads), .groups = "drop") %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = n_reads, fill = is_dmg))+
      geom_area(stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
      facet_nested(.~core*is_dmg, scales = "free_x")+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
      labs(x = "Age (y bp)", y = "Number of reads", fill = "Damage classification")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/eukaryotes/preliminary/damage_reads_species.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()

############################
## Genus level data ########
############################

gendata <- read_tsv("results/eukaryotes/taxonomy/genus_level_data.tsv")

plotA <- gendata %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg_local = fct_relevel(is_dmg_local, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = A_b, fill = is_dmg_local, size = n_reads))+
      geom_point(color = "black", shape = 21, alpha = 0.5, stroke = 0.2)+
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
      facet_nested(.~core*is_dmg)+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))+
      labs(x = "Age (y bp)", y = "Damage (A_b)", color = "Local-based classification (A_b > 0.1 & Zfit >= 2)", size = "Number of reads")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/eukaryotes/preliminary/damage_depositional_classification_genus.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()


plotA <- gendata %>%
    group_by(label, core, y_bp, is_dmg) %>%
    summarise(n_reads = sum(n_reads), .groups = "drop") %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = n_reads, fill = is_dmg))+
      geom_area(stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
      facet_nested(.~core*is_dmg, scales = "free_x")+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      # scale_y_continuous(breaks = scales::breaks_pretty(n = 3))+
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
      labs(x = "Age (y bp)", y = "Number of reads", fill = "Damage classification")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/eukaryotes/preliminary/damage_reads_genus.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()

############################
## Family level data ########
############################

famdata <- read_tsv("results/eukaryotes/taxonomy/family_level_data.tsv")

plotA <- famdata %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg_local = fct_relevel(is_dmg_local, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = A_b, fill = is_dmg_local, size = n_reads))+
      geom_point(color = "black", shape = 21, alpha = 0.5, stroke = 0.2)+
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
      facet_nested(.~core*is_dmg)+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))+
      labs(x = "Age (y bp)", y = "Damage (A_b)", color = "Local-based classification (A_b > 0.1 & Zfit >= 2)", size = "Number of reads")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/eukaryotes/preliminary/damage_depositional_classification_family.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()


plotA <- famdata %>%
    group_by(label, core, y_bp, is_dmg) %>%
    summarise(n_reads = sum(n_reads), .groups = "drop") %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = n_reads, fill = is_dmg))+
      geom_area(stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
      facet_nested(.~core*is_dmg, scales = "free_x")+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
      labs(x = "Age (y bp)", y = "Number of reads", fill = "Damage classification")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/eukaryotes/preliminary/damage_reads_family.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()




spdata %>%
  filter(is_dmg == "Damaged") %>%
  filter(species %in% c("Fratercula arctica", "Alca torda", "Emiliania huxleyi")) %>%
  select(label, species, n_reads) %>%
  complete(label, species, fill = list(n_reads = 0)) %>%
  filter(species %in% c("Fratercula arctica", "Alca torda")) %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = species))+
    geom_area(stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
    # geom_line()+
    scale_fill_manual(values = paired_genus)+
    # scale_color_manual(values = c("Fratercula arctica" = "black", "Emiliania huxleyi" = "red"))+
    # geom_smooth(method = "loess")+
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    labs(x = "Core", y = "Number of reads", fill = "Taxa")+
    facet_nested(~core, scales = "free_x")

gendata %>%
  filter(is_dmg == "Damaged") %>%
  filter(genus %in% c("Fratercula", "Alca", "Emiliania")) %>%
  select(label, genus, n_reads) %>%
  complete(label, genus, fill = list(n_reads = 0)) %>%
  filter(genus %in% c("Fratercula", "Alca")) %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = n_reads))+
    geom_area(aes(fill = genus), stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
    geom_point(data = . %>% filter(n_reads > 0), aes(fill = genus), shape = 21, color = "black")+
    scale_fill_manual(values = paired_genus)+
    scale_color_manual(values = paired_genus)+
    
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    labs(x = "Core", y = "Number of reads", fill = "Taxa")+
    facet_nested(~core, scales = "free_x")

famdata %>%
  filter(is_dmg == "Damaged") %>%
  filter(phylum %in% c("Foraminifera", "Haptophyta")) %>%
  select(label, family, n_reads) %>%
  complete(label, family, fill = list(n_reads = 0)) %>%
  filter(grepl("Glob", family)) %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = family))+
    geom_area(stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
    scale_fill_manual(values = paired_genus)+
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    labs(x = "Core", y = "Number of reads", fill = "Taxa")+
    facet_nested(~core, scales = "free_x")

famdata %>%
  filter(is_dmg == "Damaged") %>%
  filter(family %in% c("Alcidae", "Noelaerhabdaceae")) %>%
  select(label, family, n_reads) %>%
  complete(label, family, fill = list(n_reads = 0)) %>%
  filter(family %in% c("Alcidae")) %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = family))+
    geom_area(stat = "identity", position = "stack", alpha = 0.5, color = "black", linewidth = 0.3)+
    scale_fill_manual(values = paired_genus)+
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    labs(x = "Core", y = "Number of reads", fill = "Taxa")+
    facet_nested(~core, scales = "free_x")


spdata %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n_sp = n_distinct(species), .groups = "drop") %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = n_sp))+
    geom_point()+
    geom_line()+
    # geom_smooth(method = "loess")+
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    labs(x = "Core", y = "Number of species", fill = "Taxa")+
    facet_nested(~core)



spreads <- spdata %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n_taxa = n_distinct(species), .groups = "drop") %>%
  mutate(level = "species")

genreads <- gendata %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n_taxa = n_distinct(genus), .groups = "drop") %>%
  mutate(level = "genus")

famreads <- famdata %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n_taxa = n_distinct(family), .groups = "drop") %>%
  mutate(level = "family")

# rbind(spreads, genreads, famreads) %>%
#   left_join(cdata) %>%
#   mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#   ggplot(., aes(x = y_bp, y = n_reads, color = level))+
#     geom_point()+
#     geom_line()+
#     coord_flip()+
#     scale_y_sqrt()+
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
#     labs(x = "Core", y = "Number of taxa", fill = "Taxa")+
#     facet_nested(~core)

  



#### Compare taxonomic results across levels

fmA <- famdata %>%
  filter(is_dmg == "Damaged") %>%
  # filter(family == "Noelaerhabdaceae") %>% 
  # filter(label == "LV3003046952") %>% 
  group_by(core, label, family, taxa = family) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  mutate(resolution = "family")

gnA <- gendata %>%
  filter(is_dmg == "Damaged") %>%
  # filter(family == "Noelaerhabdaceae") %>%
  # filter(label == "LV3003046952") %>% 
  group_by(core, label, family, taxa = genus) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  mutate(resolution = "genus")

spA <- spdata %>% 
  filter(is_dmg == "Damaged") %>%
  # filter(family == "Noelaerhabdaceae") %>%
  # filter(label == "LV3003046952") %>% 
  group_by(core, label, family, taxa = species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  mutate(resolution = "species")
  
samples <- cdata %>% pull(label) %>% unique()

plotA <- rbind(fmA, gnA, spA) %>% View()
  group_by(label, family) %>%
  mutate(
    proportion_species = sum(n_reads[resolution == "species"]) / sum(n_reads[resolution == "family"]),
    proportion_genus = sum(n_reads[resolution == "genus"]) / sum(n_reads[resolution == "family"])) %>%
  filter(resolution == "family") %>% 
  pivot_longer(cols = c(proportion_species, proportion_genus), names_to = "taxonomic_resolution", values_to = "proportion") %>% 
  group_by(taxa) %>% 
  mutate(median_proportion = mean(proportion[taxonomic_resolution == "proportion_species"])) %>%
  ungroup() %>%
  mutate(taxa = fct_reorder(taxa, median_proportion)) %>% #select(taxa, median_proportion) %>% distinct() %>% View()
  left_join(famdata %>% select(kingdom, phylum, class, order, family) %>% distinct()) %>%
  mutate(taxonomic_resolution = fct_relevel(taxonomic_resolution, "proportion_species", "proportion_genus"))  %>%
  ggplot(., aes(x = taxa, y = proportion, fill = paste(kingdom, phylum)))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(aes(size = n_reads))+
    scale_size_continuous(range = c(0.5, 3))+
    coord_flip()+
    scale_fill_manual(values = paired_genus)+
    facet_nested(~taxonomic_resolution)+
    # facet_nested(kingdom*phylum~., scales = "free_y", space = "free_y")+
    guides(fill = guide_legend(ncol = 1))+
    labs(x = "Taxa", y = "Proportion of reads assigned at species level", fill = "Taxa")
    # theme(legend.position = "bottom")
plotA



rbind(fmA, gnA, spA) %>%
  select(-taxa) %>% 
  filter(family == "Noelaerhabdaceae") %>%
  group_by(core, label, family, resolution) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  rename(target_family = family) %>% 
  pivot_wider(names_from = resolution, values_from = n_reads, values_fill = 0) %>%
  mutate(proportion_species = species / family) %>% 
  left_join(cdata) %>% View()
  ggplot(., aes(x = y_bp))+
    geom_point(aes(y = family), color = "blue", size = 3)+
    geom_point(aes(y = species), color = "red", size = 3)+
    geom_segment(aes(y = family, yend = species))+
    coord_flip()+
    facet_nested(~core, scales = "free_x")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))
  
  
  
  
  mutate(
    proportion_species = sum(n_reads[resolution == "species"]) / sum(n_reads[resolution == "family"]),
    proportion_genus = sum(n_reads[resolution == "genus"]) / sum(n_reads[resolution == "family"])) %>%
  filter(resolution == "family") %>% 
  pivot_longer(cols = c(proportion_species, proportion_genus), names_to = "taxonomic_resolution", values_to = "proportion") %>% 
  group_by(taxa) %>% 
  mutate(median_proportion = mean(proportion[taxonomic_resolution == "proportion_species"])) %>%
  ungroup() %>%
  mutate(taxa = fct_reorder(taxa, median_proportion)) %>% #select(taxa, median_proportion) %>% distinct() %>% View()
  left_join(famdata %>% select(kingdom, phylum, class, order, family) %>% distinct()) %>%
  mutate(taxonomic_resolution = fct_relevel(taxonomic_resolution, "proportion_species", "proportion_genus"))  %>% head()
  ggplot(., aes(x = taxa, y = proportion, fill = paste(kingdom, phylum, class)))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(aes(size = n_reads))+
    scale_size_continuous(range = c(0.5, 3))+
    coord_flip()+
    scale_fill_manual(values = paired_genus)+
    facet_nested(~taxonomic_resolution)+
    # facet_nested(kingdom*phylum~., scales = "free_y", space = "free_y")+
    guides(fill = guide_legend(ncol = 1))+
    labs(x = "Taxa", y = "Proportion of reads assigned at species level", fill = "Taxa")
    # theme(legend.position = "bottom")
plotA




rbind(fmA, gnA, spA) %>%
  group_by(label, family) %>%
  mutate(
    proportion_species = sum(n_reads[resolution == "species"]) / sum(n_reads[resolution == "family"]),
    proportion_genus = sum(n_reads[resolution == "genus"]) / sum(n_reads[resolution == "family"])) %>%
  filter(resolution == "family") %>% 
  pivot_longer(cols = c(proportion_species, proportion_genus), names_to = "taxonomic_resolution", values_to = "proportion") %>% 
  group_by(taxa) %>% 
  mutate(median_proportion = mean(proportion[taxonomic_resolution == "proportion_species"])) %>%
  ungroup() %>%
  mutate(taxa = fct_reorder(taxa, median_proportion)) %>% #select(taxa, median_proportion) %>% distinct() %>% View()
  left_join(famdata %>% select(kingdom, phylum, class, order, family) %>% distinct()) %>%
  mutate(taxonomic_resolution = fct_relevel(taxonomic_resolution, "proportion_species", "proportion_genus")) %>%
  left_join(cdata) %>% 
  filter(median_proportion > 0.3 & median_proportion < 0.8) %>%
  ggplot(., aes(x = y_bp, y = proportion, color = paste(kingdom, phylum, class)))+
    geom_point()+
    geom_line()+
    # geom_smooth(method = "loess")+
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    labs(x = "Age (y bp)", y = "Proportion of reads assigned at species level", color = "Taxa")+
    theme(legend.position = "bottom")+
    facet_nested(~taxonomic_resolution*core)









Rspdata %>%
  filter(family_unicorn == "Noelaerhabdaceae", library_id == "LV3003046952") %>%
  summarise(
    total_reads_from_unicorn = sum(total_reads_from_unicorn), 
    total_n_reads_from_mD = sum(total_n_reads_from_mD), 
    .groups = "drop")

Rgendata %>% 
  filter(family_unicorn == "Noelaerhabdaceae", library_id == "LV3003046952") %>%
  summarise(
    total_reads_from_unicorn = sum(total_reads_from_unicorn), 
    total_n_reads_from_mD = sum(total_n_reads_from_mD), 
    .groups = "drop")

Rfamdata %>% filter(library_id == "LV7008867186") %>% View()
  summarise(
    total_reads_from_unicorn = sum(total_reads_from_unicorn), 
    total_n_reads_from_mD = sum(total_n_reads_from_mD), 
    .groups = "drop")
Rfamdata %>% 
    filter(family.x_unicorn == "Noelaerhabdaceae", library_id == "LV3003046952") %>%
    summarise(
      total_reads_from_unicorn = sum(total_reads_from_unicorn), 
      total_n_reads_from_mD = sum(total_n_reads_from_mD), 
      .groups = "drop")

spdata %>%
  filter(family == "Noelaerhabdaceae", label == "LV3003046952") %>%
  summarise(total_n_reads = sum(n_reads), .groups = "drop")

famdata %>%
  filter(family == "Noelaerhabdaceae", label == "LV3003046952") %>%
  summarise(total_n_reads = sum(n_reads), .groups = "drop")

  head()
