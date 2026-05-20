library(gghighlight)
library(tidyverse)
library(biomformat)
library(phyloseq)
library(showtext)
library(ggh4x)
library(janitor)
library(ggvenn)
library(dendextend)
library(scales)

setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")
source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")

cdata <- read_tsv(file = "data/metadata_v5.tsv")
tax_data <- read_tsv("results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz")
reference_summary <- read_tsv("results/microbial/damage/damage-classification-depositional/reference-summary.tsv")
sample_baselines <- read_tsv("results/microbial/damage/damage-classification-depositional/sample-baselines.tsv")


# smp <- read_tsv("/projects/caeg/people/ngm902/apps/repos/aeDNA/config/samples_pp.tsv")
# units <- read_tsv("/projects/caeg/people/ngm902/apps/repos/aeDNA/config/units_pp.tsv")

# cdata %>%
#   filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>% nrow()

# smp %>%
#   filter(sample %in% cdata[cdata$core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"), "label"]$label) %>%
#   write_tsv("/projects/caeg/people/ngm902/apps/repos/aeDNA/config/samples_rocs.tsv")
# units %>%
#   filter(sample %in% cdata[cdata$core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"), "label"]$label) %>%
#   write_tsv("/projects/caeg/people/ngm902/apps/repos/aeDNA/config/units_rocs.tsv")


tax_data %>% 
  group_by(is_dmg) %>%
  count() %>%
  ungroup() %>%
  rename(number = n) %>%
  mutate(percent = number/sum(number))

plotA <- tax_data %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg_local = fct_relevel(is_dmg_local, "Non-damaged",  "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = A_b, fill = is_dmg_local, size = n_reads))+
      geom_point(color = "black", shape = 21, alpha = 0.5, stroke = 0.2)+
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
      facet_nested(.~core*is_dmg)+
      scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
      # scale_y_continuous(breaks = scales::breaks_pretty(n = 3))+
      scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))+
      labs(x = "Age (y bp)", y = "Damage (A_b)", color = "Local-based classification (A_b > 0.1 & Zfit >= 2)", size = "Number of reads")+
      coord_flip()+
      theme(legend.position = "bottom")

png(file = "plots/microbial/preliminary/damage_depositional_classification.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()


plotA <- tax_data %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
    ggplot(., aes(x = A_b, y = y_bp, color = is_dmg))+
        geom_point(size = 1, alpha = 0.6)+
        geom_point(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), sd_A_b = sd(A_b), .groups = "drop"), 
            aes(x = mean_A_b, fill = is_dmg), size = 3, shape = 21, color = "black") +
        geom_errorbarh(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), n = sum(!is.na(A_b)), sd_A_b = sd(A_b)/sqrt(n), .groups = "drop"), 
            aes(x = mean_A_b, xmin = mean_A_b - sd_A_b, xmax = mean_A_b + sd_A_b), color = "black", height = 5, alpha = 1) +
        labs(x = "Damage (A_b)", y = "Age (y bp)", fill = "References classified as", color = "References classified as")+
        scale_y_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
        facet_nested(~core) +
        theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

png(file = "plots/microbial/preliminary/damage_time_median.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()




plotA <- tax_data %>%
    group_by(label, is_dmg, domain) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    complete(label, is_dmg, domain, fill = list(abundance = 0)) %>%
    left_join(cdata) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
    ggplot(., aes(x = y_bp, y = abundance, fill = domain))+
        geom_area(position = "stack", stat = "identity", color = "black", linewidth = 0.2)+
        facet_nested(.~core*is_dmg, scales = "free_x")+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::breaks_pretty(n = 3))+
        scale_fill_manual(values= c("d__Viruses"="#f3c96d", "d__Bacteria"="#74afb9", "d__Archaea"="#9f443d"))+
        guides(fill = guide_legend(nrow = 1))+
        labs(x = "Age (y bp)", y = "Abundance", fill = "Domain")+
        coord_flip()+
        theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

png(file = "plots/microbial/preliminary/abundance_domain.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()







to_plot <- tax_data %>%
  mutate(tax = paste(domain, phylum, sep = "; ")) %>%
  group_by(label, is_dmg, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(label, is_dmg, tax, fill = list(abundance = 0)) %>%
  group_by(label, is_dmg) %>%
  mutate(abundance = abundance/sum(abundance)) %>%
  ungroup() %>%

  group_by(tax) %>%
  mutate(total_abundance = mean(abundance)) %>%
  ungroup() %>%

  mutate(tax = ifelse(total_abundance >= 0.005, tax, "Other (rel. abund < 0.005)")) %>%

  group_by(label, is_dmg, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(label, is_dmg, tax, fill = list(abundance = 0)) %>%

  group_by(tax) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%

  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
  mutate(tax = {
    tmp <- fct_reorder(tax, -total_abundance)
    other_lvls <- levels(tmp)[grepl("^Other", levels(tmp))]
    if (length(other_lvls) > 0) fct_relevel(tmp, other_lvls, after = Inf) else tmp
  })

cc <- paired_genus[1:length(levels(to_plot$tax))]
names(cc) <- levels(to_plot$tax)
cc[grepl("^Other", names(cc))] <- "grey"

plotA <- to_plot %>%
  ggplot(., aes(x = y_bp, y = abundance, fill = tax))+
    geom_area(position = "stack", stat = "identity", color = "black", linewidth = 0.2)+
    facet_nested(.~core*is_dmg, scales = "free_x")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::breaks_pretty(n = 3))+
    scale_fill_manual(values = cc)+
    guides(fill = guide_legend(nrow = 4, title = NULL))+
    labs(x = "Age (y bp)", y = "Relative abundance", size = "Number of reads")+
    coord_flip()+
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

plotB <- tax_data %>%
  mutate(tax = paste(domain, phylum, sep = "; ")) %>%
  group_by(is_dmg, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(is_dmg, tax, fill = list(abundance = 0)) %>%
  group_by(is_dmg) %>%
  mutate(abundance = abundance/sum(abundance)) %>%
  ungroup() %>%

  group_by(tax) %>%
  mutate(mean_abundance = mean(abundance)) %>%
  ungroup() %>%

  pivot_wider(names_from = is_dmg, values_from = abundance) %>%
  janitor::clean_names() %>%
  mutate(ratio = log2((damaged + 1e-6)/(non_damaged + 1e-6))) %>%

  mutate(to_color = ifelse(tax %in% to_plot$tax, tax, "Other (rel. abund < 0.005)")) %>%
  mutate(phylum = str_extract(tax, "p__[^_]+")) %>%
  ggplot(., aes(x = ratio, y = mean_abundance, fill = to_color))+
    geom_point(shape = 21, color = "black", size = 4)+
    ggrepel::geom_text_repel(aes(label = ifelse(mean_abundance > 0.02 | ratio > 3 | ratio < -5 | !grepl("Other", to_color), phylum, "")), size = 3, color = "grey28")+
    
    scale_fill_manual(values = cc)+
    guides(fill = guide_legend(nrow = 2))+
    
    labs(x = expression(log[2](Damaged / Non-damaged~relative~abundance)), y = "Mean total relative abundance", fill = "Taxon", title = "")+
    
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
    theme(legend.position = "none", aspect.ratio = 1)
  
png(file = "plots/microbial/preliminary/phyla_damage_status.png", width = 15.3, height = 20, units = "in", res = 600)
egg::ggarrange(plotA, plotB, ncol = 1, heights = c(2.6, 1))
dev.off()







plotA <- tax_data %>%
  group_by(label, is_dmg) %>%
  summarise(
    # subspecies = n_distinct(subspecies),
    species = n_distinct(species),
    genus = n_distinct(genus),
    family = n_distinct(family),
    # n_order = n_distinct(order),
    # n_class = n_distinct(class),
    # n_phylum = n_distinct(phylum),
    .groups = "drop") %>%
  pivot_longer(cols = c(species, genus, family), names_to = "taxonomic_level", values_to = "n_taxa") %>%
  left_join(cdata) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
  ggplot(., aes(x = y_bp, y = n_taxa, color = taxonomic_level))+
    geom_point(size = 1)+
    geom_line(linewidth = 0.5)+
    labs(x = "Age (y bp)", y = "Number of taxa", color = "Taxonomic level")+
    coord_flip()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
    facet_nested(~core*is_dmg) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

png(file = "plots/microbial/preliminary/n_taxa_damage_status.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()








tax_data %>%
  mutate(tax = paste(domain, phylum)) %>%
  group_by(is_dmg, core, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  group_by(core, is_dmg) %>%
  mutate(abundance = abundance/sum(abundance)) %>%
  ungroup() %>%

  group_by(tax) %>%
  mutate(total_abundance = mean(abundance)) %>%
  ungroup() %>%

  mutate(tax = ifelse(total_abundance >= 0.001, tax, "Other")) %>%

  group_by(core, is_dmg, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(core, is_dmg, tax, fill = list(abundance = 0)) %>%

  group_by(tax) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%

  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
  ggplot(., aes(x = core, y = abundance, fill = tax))+
    geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.2)+
    facet_nested(.~is_dmg, scales = "free_x")+
    coord_flip()+
    guides(fill = guide_legend(ncol = 1))+
    scale_fill_manual(values = paired_genus)+
    labs(x = "Age (y bp)", y = "Abundance", fill = "Domain", size = "Number of reads")





tax_data %>%
  mutate(tax = paste(domain, phylum)) %>%
  group_by(is_dmg, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  group_by(is_dmg) %>%
  mutate(abundance = abundance/sum(abundance)) %>%
  ungroup() %>%
  complete(is_dmg, tax, fill = list(abundance = 0)) %>%
  pivot_wider(names_from = is_dmg, values_from = abundance) %>%
  janitor::clean_names() %>%
  mutate(
    mean_abundance = (damaged + non_damaged)/2,
    ratio = log2((damaged + 1e-6)/(non_damaged + 1e-6))
  ) %>%
  mutate(tax = fct_reorder(tax, ratio)) %>%
  filter(mean_abundance > 0.0001) %>%
  ggplot(., aes(x = tax, y = ratio, size = mean_abundance))+
    geom_point()+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
     scale_size_continuous(range = c(1, 10), breaks = c(1e-6, 1e-4, 1e-2, 1))+
     scale_y_continuous(breaks = scales::breaks_pretty(n = 5))+
    #  scale_x_discrete(labels = function(x) gsub("d__|p__", "", x))+
    labs(x = "Taxon (Domain_Phylum)", y = "Log2 ratio of abundance (Damaged/Non-damaged)", fill = "Taxon")+
    coord_flip()




  group_by(tax) %>%
  mutate(total_abundance = mean(abundance)) %>%
  ungroup() %>%

  mutate(tax = ifelse(total_abundance >= 0.001, tax, "Other")) %>%

  group_by(core, is_dmg, tax) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%

  group_by(tax) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%

  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged",  "Damaged")) %>%
  ggplot(., aes(x = core, y = abundance, fill = tax))+
    geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.2)+
    facet_nested(.~is_dmg, scales = "free_x")+
    coord_flip()+
    guides(fill = guide_legend(ncol = 1))+
    scale_fill_manual(values= c("d__Viruses"="#f3c96d", "d__Bacteria"="#74afb9", "d__Archaea"="#9f443d"))+
    scale_fill_manual(values = paired_genus)+
    labs(x = "Age (y bp)", y = "Abundance", fill = "Domain", size = "Number of reads")

