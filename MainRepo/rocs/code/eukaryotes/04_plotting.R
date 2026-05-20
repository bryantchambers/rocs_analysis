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
cdata <- read_tsv(file = "results/common/metadata_xrf_matched.tsv")

# ------------------------------
# Core analysis settings
# ------------------------------
core_levels <- c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")
taxonomy_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

species <- readRDS("results/eukaryotes/taxonomy/species_exploration.rds")
genus <- readRDS("results/eukaryotes/taxonomy/genus_exploration.rds")
family <- readRDS("results/eukaryotes/taxonomy/family_exploration.rds")




species$summaries$sample_summary_dmg %>%
  left_join(cdata %>% select(-y_bp, -depth_in_core_cm, -core), by = "label") %>%
  ggplot(., aes(x = y_bp, y = observed_richness)) +
    geom_point() +
    coord_flip() +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    facet_nested(~core)

species$summaries$sample_summary_dmg %>%
  left_join(cdata %>% select(-y_bp, -depth_in_core_cm, -core), by = "label") %>%
  ggplot(., aes(x = y_bp, y = Shannon)) +
    geom_point() +
    coord_flip() +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    facet_nested(~core)




genus$summaries$group_summary_dmg %>%
  select(label, ecological_group, relative_abundance) %>%
  complete(label, ecological_group, fill = list(relative_abundance = 0)) %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = relative_abundance, fill = ecological_group)) +
    geom_area(stat = "identity", position = "stack") +
    coord_flip() +
    scale_fill_manual(values = paired_genus)+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    facet_nested(~core)




species$data$dmg %>% 
  filter(ecological_group == "Marine mammals") %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = n_reads, fill = species)) +
      geom_area(stat = "identity", position = "stack") +
      coord_flip() +
      scale_fill_manual(values = paired_genus)+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      facet_nested(~core)
