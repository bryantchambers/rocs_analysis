library(lvplot)
library(showtext)
library(tidyverse)
library(DescTools)
library(ggh4x)
library(gghighlight)
library(data.table)
library(dtplyr)

# showtext_auto()
# httpgd::hgd()
# httpgd::hgd_browse()


# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/get-metadata.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/dmg.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/get_calculate_plot_grid.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/get_penalized_weighted_median_reads.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/perk.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/damage_est_function.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/perk_wrapper.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/perk_wrapper_function.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/mikkel_euk/get_dmg_decay_fit.R")

source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")
source("/projects/caeg/people/ngm902/rocs/v3/code/lib_rocs.R")

setwd("/projects/caeg/people/ngm902/rocs/v3")


# Let's load the cdata
cdata <- read.table(file = "associated_data/metadata_v3.txt", sep = "\t", header = T)


## Load data  ##
subspecies_stats_mapping <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")





plotA <- subspecies_stats_mapping %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(n = n_distinct(subspecies), n_reads = sum(n_reads), .groups = "drop") %>%
  left_join(cdata) %>%
  mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
  filter(core == "GeoB25202_R1" | core == "GeoB25202_R2") %>%
  ggplot(., aes(x = y_bp, y = n_reads))+
    geom_point()+
    geom_line()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(x = "Age (y bp)", y = "Number of reads", title = "(Prok) Number of reads")+
    facet_nested(~core, scales = "free")

plotB <- subspecies_stats_mapping %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(n = n_distinct(subspecies), n_reads = sum(n_reads), .groups = "drop") %>%
  left_join(cdata) %>%
  mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
  filter(core == "GeoB25202_R1" | core == "GeoB25202_R2") %>%
  ggplot(., aes(x = y_bp, y = n))+
    geom_point()+
    geom_line()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))+
    coord_flip()+
    labs(x = "", y = "Number of species", title = "Number of species")+
    facet_nested(~core, scales = "free")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Calculate Shannon diversity index per sample (label)
shannon_diversity <- subspecies_stats_mapping %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label) %>%
  summarise(
    shannon = -sum(ifelse(n_reads > 0, (n_reads / sum(n_reads)) * log(n_reads / sum(n_reads)), 0)),
    n = n_distinct(subspecies),
    n_reads = sum(n_reads),
    .groups = "drop"
  ) %>%
  left_join(cdata) %>%
  mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
  filter(core == "GeoB25202_R1" | core == "GeoB25202_R2")

# Plot Shannon diversity
plotC <- ggplot(shannon_diversity, aes(x = y_bp, y = shannon)) +
  geom_point() +
  geom_line() +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4)) +
  coord_flip() +
  labs(x = "", y = "Shannon diversity", title = "Shannon diversity") +
  facet_nested(~core, scales = "free") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plotD <- subspecies_stats_mapping %>%
  filter(is_dmg == "Damaged") %>%
  group_by(label, phylum) %>%
  summarise(n_reads = sum(abundance), .groups = "drop") %>%
  complete(label, phylum, fill = list(n_reads = 0)) %>%
  group_by(label) %>%
  mutate(n_reads = n_reads / sum(n_reads)) %>%
  ungroup() %>%
  filter(phylum %in% c("p__Hadarchaeota", "p__Aerophobota")) %>%
  # filter(phylum %in% c("p__Asgardarchaeota", "p__Hadarchaeota", "p__Aerophobota", "p__Atribacterota", "p__KSB1", "p__UBA6262", "p__UBP18")) %>%
  left_join(cdata) %>%
  mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
  filter(core == "GeoB25202_R1" | core == "GeoB25202_R2") %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = phylum))+
    geom_area(stat = "identity", position = "stack") +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(x = "", y = "Proportion of abundance", title = "Deep sediment taxa")+
    facet_nested(~core, scales = "free")+
    scale_fill_manual(values = paired_genus)+
    guides(fill = guide_legend(ncol = 1)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "/projects/caeg/people/ngm902/rocs/v3/results/plots/microbial/preliminary_exp_prok.png", width = 17, height = 10, units = "in", res = 600)
egg::ggarrange(plotA, plotB, plotC, plotD, nrow = 1)
dev.off()
















filtered_subspecies <- subspecies_stats_mapping %>%
  group_by(phylum, class, order, family, genus, species, subspecies, is_dmg) %>%
  summarise(n_sites_100 = n(), tota_reads_100 = sum(n_reads), median_reads_100 = median(n_reads), .groups = "drop") %>%
  group_by(phylum, class, order, family, genus, species, subspecies) %>%
  mutate(
    total_sites = sum(n_sites_100),
    proportion_good_fit = if (any(is_dmg == "Damaged")) n_sites_100[is_dmg == "Damaged"] / total_sites else 0
  ) %>% 
  ungroup() %>%
  filter(proportion_good_fit != 0)

subspecies_stats_mapping <- subspecies_stats_mapping %>%
  mutate(selected = case_when(
    is_dmg == "Damaged" & subspecies %in% as.character(filtered_subspecies %>% filter(is_dmg == "Damaged" & n_sites_100 >= 5) %>% pull(subspecies)) ~ "yes",  # Damaged observations with an occurrence threshold of >= 5
    subspecies %in% as.character(filtered_subspecies %>% filter(is_dmg == "Damaged" & n_sites_100 >= 10 & proportion_good_fit > 0.5) %>% pull(subspecies)) ~ "yes", # All species with >0.5 proportion of damaged
    TRUE ~ "no"
  ))

subspecies_stats_mapping %>% 
  write_tsv("results/microbial/taxonomy/dmg-summary-ssp_selected.tsv.gz")


subspecies_stats_mapping %>%
  group_by(label, selected, is_dmg) %>%
  summarise(
    abundance = sum(abundance), 
    n_reads = sum(n_reads), 
    .groups = "drop"
  ) %>%
  complete(label, selected, is_dmg, fill = list(abundance = 0, n_reads = 0)) %>%
  left_join(cdata) %>%
  ggplot(., aes(y = abundance, x = y_bp, fill = is_dmg))+
    geom_area(stat = "identity", position = "stack")+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e5))+
        scale_y_sqrt()+
        facet_nested(~core*selected, scales = "free") +
        coord_flip()


filtered_subspecies %>%
  filter(!subspecies %in% aa$subspecies) %>%
  # filter(fit == "bad") %>% 
  View()

filt <- filtered_subspecies %>% 
      filter(fit == "good", n_sites_100 >= 5, proportion_good_fit > 0.5)

subspecies_stats_mapping %>%
      mutate(kept = ifelse(subspecies %in% filt$subspecies, "yes", "no")) %>%
      group_by(kept, label, y_bp, core) %>%
      summarise(
        abundance = sum(abundance), 
        n_reads = sum(n_reads), 
        .groups = "drop"
      ) %>% 
      complete(label, kept, fill = list(abundance = 0, n_reads = 0)) %>%
      left_join(cdata) %>%
      ggplot(., aes(x = y_bp, y = abundance, fill = kept)) +
        geom_area(stat = "identity", position = "stack")+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e3))+
        facet_nested(~core, scales = "fixed") +
        coord_flip()

subspecies_stats_mapping %>%
      mutate(kept = ifelse(subspecies %in% filt$subspecies, "yes", "no")) %>%
      group_by(kept, label, y_bp, core) %>%
      summarise(
        abundance = sum(abundance), 
        n_reads = sum(n_reads), 
        .groups = "drop"
      ) %>%
      group_by(kept, core) %>%
      summarise(abundance = sum(abundance), .groups = "drop") %>%
      group_by(core) %>%
      mutate(rel_abundance = abundance / sum(abundance)) %>%
      ungroup() %>%
      print(n = 30)

subspecies_stats_mapping %>%
      group_by(phylum, class, order, family, genus, species, subspecies) %>%
      summarise(
        n = n_distinct(label), 
        abundance = sum(abundance), 
        n_reads = sum(n_reads), 
        .groups = "drop"
      ) %>% 
      mutate(kept = ifelse(subspecies %in% filt$subspecies, "yes", "no")) %>%
      mutate(phylum = paste(phylum, class, sep = "; ")) %>%
      group_by(phylum, kept) %>%
      summarise(abundance = sum(abundance), .groups = "drop") %>%
      group_by(kept) %>%
      mutate(rel_abundance = abundance / sum(abundance)) %>%
      ungroup() %>%
      group_by(phylum) %>%
      mutate(total = sum(rel_abundance)) %>%
      ungroup() %>%
      mutate(phylum = fct_reorder(phylum, -total)) %>%
      ggplot(., aes(x = kept, y = rel_abundance, fill = phylum)) +
        geom_bar(stat = "identity", position = "stack")+
        scale_fill_manual(values = paired_genus)+
        theme(legend.position = "bottom")+
        coord_flip()

filtered_subspecies %>% pull(subspecies) %>% unique() %>% length()
fs <- filtered_subspecies %>%
  filter(fit == "good", n_sites_100 >= 5, proportion_good_fit > 0.5) %>%
  pull(subspecies) #%>% unique() %>% length()


filtered_subspecies %>%
  ggplot(., aes(x = n_sites_100, y = proportion_good_fit, size = tota_reads_100, color = ifelse(subspecies %in% fs, "yes", "no"))) +
  geom_point()

filtered_subspecies %>%
  filter(fit == "good" & n_sites_100 >= 5 & proportion_good_fit > 0.15) %>%
  pull(species)

subspecies_stats_mapping %>%
    mutate(is_dmg = ifelse(y_bp < 10000 & fit == "good", "Damaged", is_dmg)) %>%
    filter(!grepl("GeoB25206|ST5", core)) %>%
    filter(is_dmg == "Damaged") %>%
    select(label, domain, phylum, class, order, family, genus, species, subspecies, abundance, is_dmg, A_b) %>%
    left_join(cdata)


taxatable <- subspecies_stats_mapping %>%
  select(domain, phylum, class, order, family, genus, species, subspecies) %>%
  distinct() %>%
  arrange(domain, phylum, class, order, family, genus, species, subspecies)

subspecies_stats <- subspecies_stats_mapping %>%
  group_by(domain, phylum, class, order, family, genus, species, subspecies) %>%
  summarise(abundance = sum(abundance), n = n(), .groups = "drop") %>%
  arrange(domain, phylum, class, order, family, genus, species, subspecies)


cdata %>% filter(y_bp < 100000 & y_bp > 90000, core == "ST13")
subset <- subspecies_stats_mapping %>%
  filter(
    y_bp < 120000 & y_bp > 80000,
    core == "ST13",
    is_dmg == "Damaged") 

plotB <- subset %>% 
  group_by(label) %>%
  summarise(count = n(), .groups = "drop")

plotB <- subset %>% 
  group_by(label, species) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(species = unique(species),
        label = unique(label),
        fill = list(abundance = 0)) %>%
  left_join(cdata) %>%
  left_join(subspecies_stats_mapping %>% select(domain, phylum, class, order, family, genus, species) %>% distinct()) %>%
  group_by(label, core, y_bp, family) %>%
  summarise(abundance = sum(abundance), n = n(), .groups = "drop") %>%
  # mutate(color = ifelse(label %in% c("LV3003061474", "LV7008886768", "LV3003061477", "LV7008886715", "LV7008867130"), TRUE, FALSE)) %>%
  ggplot(., aes(x = y_bp, y = abundance, fill = family))+
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 200) +
    # geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    # geom_vline(data = . %>% filter(color == TRUE), aes(xintercept = y_bp), linetype = "dashed", color = "black") +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e3))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    scale_fill_manual(values = paired_genus) +
    labs(fill = "")
plotB



subspecies_stats_mapping %>%
  group_by(label, core, y_bp, is_dmg) %>%
  summarise(abundance = sum(abundance), n = n(), .groups = "drop") %>%
  filter(y_bp < 140000, is_dmg == "Damaged", core == "ST13") %>%
  ggplot(., aes(x = y_bp, y = abundance))+
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")+
    theme(aspect.ratio = 3/1)

## Additional plots
plottemp <- subspecies_stats_mapping %>%
	select(y_bp, temp) %>%
	distinct() %>%
	ggplot(., aes(x = y_bp, y = temp))+
		geom_line()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
		coord_flip()+
		labs(x = "Age (y bp)", y = "SST (ºC)", title = "SST")

plotmi <- subspecies_stats_mapping %>%
	select(y_bp, mis) %>%
	distinct() %>%
	ggplot(., aes(x = y_bp, y = mis))+
		geom_line()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
		scale_y_reverse()+
		coord_flip()+
		labs(x = "", y = "Benthic d18O", title = "MIS")+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

totalReads <- subspecies_stats_mapping %>%
  group_by(label) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  left_join(cdata) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(aes(x = y_bp, y = abundance)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")+
    theme(aspect.ratio = 3/1)



### Preliminary view of taxa
taxa <- subspecies_stats_mapping %>%
  filter(
    core == "ST13",
    y_bp < 120000,
  ) %>%
  group_by(species) %>%
  summarise(abundance = sum(abundance), n = n(), .groups = "drop") %>%
  arrange(desc(n))

datafilt <- subspecies_stats_mapping %>%
  # Filter data
  filter(
    y_bp < 120000,
    core == "ST13",
    species %in% taxa$species[taxa$abundance >= 1000 & taxa$n >= 5],
    # label %in% samples$label[samples$n_reads >= 500 & samples$n >= 3],
  ) %>% 
  # Complete the data
  group_by(label, species) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(abundance = 0)) %>%
  # Add taxa and core information
  left_join(subspecies_stats_mapping %>% select(domain, phylum, class, order, family, genus, species) %>% distinct()) %>%
  left_join(cdata)

plotA <- datafilt %>%  
  group_by(species) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  group_by(phylum) %>%
  mutate(phylum_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -total_abundance)) %>%
  mutate(phylum = fct_reorder(phylum, -phylum_abundance)) %>%
  ggplot(aes(x = y_bp, y = abundance, fill = phylum))+
			geom_area(stat = "identity", position = "stack", alpha = 0.9, color = "black", linewidth = 0.1) +
			# labs(title = title, x = "", y = "Number of reads", fill = "") +
      facet_nested(~domain, scales = "free") +
			# facet_nested(~core, scales = "free") +
			guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
			# scale_fill_manual(values = paired_genus) +
			scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
			scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
			coord_flip()+
			theme(aspect.ratio = 4/1, legend.position = "none")

png(file = "results/plots/test.png", width = 50, height = 10, units = "in", res = 600)
plotA
dev.off()

datafilt %>% filter(grepl("Colwellia", species))

plotA <- datafilt %>%  
  filter(phylum == "p__Nitrospirota") %>%
  # filter(order %in% c("o__Pelagibacterales", "o__Pseudomonadales", "o__Enterobacterales", "o__Flavobacteriales") ) %>%
  # filter(genus %in% c("g__Arctic96AD-7", "g__Polaribacter", "g__Colwellia_A", "g__Colwellia_A") ) %>%
  group_by(species) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(phylum_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -total_abundance)) %>%
  mutate(class = fct_reorder(class, -phylum_abundance)) %>%
  ggplot(aes(x = y_bp, y = abundance, fill = species))+
			geom_area(stat = "identity", position = "stack", alpha = 0.9, color = "black", linewidth = 0.1) +
			# labs(title = title, x = "", y = "Number of reads", fill = "") +
      facet_nested(~species, scales = "free") +
			# facet_nested(~core, scales = "free") +
			guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
			scale_fill_manual(values = paired_genus) +
			scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
			scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
			coord_flip()

png(file = "results/plots/test.png", width = 30, height = 10, units = "in", res = 600)
plotA
dev.off()



taxatable %>%
  filter(phylum == "p__Methylomirabilota")


subspecies_stats_mapping %>%
  filter(core == "GeoB25202_R1") %>%
  filter(y_bp > 120000 & y_bp < 150000) %>%
  filter(is_dmg == "Damaged") %>%
  head()






# Create consistent color palettes
taxa_all <- subspecies_stats_mapping %>%
  filter(
    core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"),
    # y_bp > x_min_s-25000 & y_bp < x_max_s+20000
  ) %>%
  group_by(species) %>%
  summarise(abundance = sum(abundance), n = n(), .groups = "drop") %>%
  filter(
    abundance >= 1000 & n >= 5
  ) %>%
  arrange(desc(abundance)) %>%
  left_join(subspecies_stats_mapping %>% select(domain, phylum, class, order, family, genus, species) %>% distinct())

# 1
names_tmp <- taxa_all %>% filter(phylum %in% c("p__Cyanobacteria")) %>% pull(species)
colors_Bacillariophyta <- paired_genus[1:length(names_tmp)]
names(colors_Bacillariophyta) <- names_tmp

# 2
names_tmp <- taxa_all %>% filter(phylum %in% c("p__Methylomirabilota")) %>% pull(species)
colors_Chlorophyta <- paired_genus[1:length(names_tmp)]
names(colors_Chlorophyta) <- names_tmp

# 3
names_tmp <- taxa_all %>% filter(phylum == "Haptophyta") %>% pull(species)
colors_Haptophyta <- paired_genus[1:length(names_tmp)]
names(colors_Haptophyta) <- names_tmp

# 4
names_tmp <- taxa_all %>% filter(class %in% c("Pelagophyceae", "Bolidophyceae", "Cryptophyceae", "Dictyochophyceae", "Eustigmatophyceae", "Polycystinea")) %>% pull(species)
colors_opp <- paired_genus[1:length(names_tmp)]
names(colors_opp) <- names_tmp

# 5
names_tmp <- taxa_all %>% filter(phylum %in% c("Arthropoda", "Foraminifera", "Cercozoa", "Priapulida")) %>% pull(species)
colors_arthopoda <- paired_genus[1:length(names_tmp)]
names(colors_arthopoda) <- names_tmp

# 6
names_tmp <- taxa_all %>% filter(class %in% c("Actinopteri", "Chondrichthyes")) %>% pull(species)
colors_fishes <- paired_genus[1:length(names_tmp)]
names(colors_fishes) <- names_tmp

# 7
names_tmp <- taxa_all %>% filter(class %in% c("Mammalia")) %>% pull(species)
colors_mammalia <- paired_genus[1:length(names_tmp)]
names(colors_mammalia) <- names_tmp






# Set parameters for the plots
x_max_s <- 120000
x_min_s <- 0
core_target  <- c("ST13")
trans_applied <- "sqrt"
min_n <- 3
min_reads <- 1000


# Filter table based on samples and taxa
samples <- species_stats_mapping %>%
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    core %in% core_target
  ) %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))

taxa <- species_stats_mapping %>%
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    core %in% core_target
  ) %>%
  group_by(species) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))



datafilt <- species_stats_mapping %>%
  # Filter data
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+25000,
    species %in% taxa$species[taxa$n_reads >= min_reads & taxa$n >= min_n],
    # label %in% samples$label[samples$n_reads >= 500 & samples$n >= 3],
    core %in% core_target
  ) %>% 
  # Complete the data
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  # Add taxa and core information
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  left_join(cdata)
 

## Primary producers
plotBacillariophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Bacillariophyta",
  trans = trans_applied, y_labels = TRUE, 
  core_levels = core_target, title = "Bacillariophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Bacillariophyta)

plotChlorophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Chlorophyta",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Chlorophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Chlorophyta)

plotHaptophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Haptophyta",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Haptophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Haptophyta)

plotOthers <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = c("Pelagophyceae", "Bolidophyceae", "Cryptophyceae", "Dictyochophyceae", "Eustigmatophyceae", "Polycystinea"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Other primary producers",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_opp)

## Intermediary and top consumers
plotArthropoda <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = c("Arthropoda", "Foraminifera", "Cercozoa", "Priapulida"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Microorg. consummers",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_arthopoda)

plotActinopteri <- plot_reads_by_taxa(datafilt,
  taxa_col = "class", taxa_value = c("Actinopteri", "Chondrichthyes"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Fishes",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_fishes)

plotMammalia <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Mammalia",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Mammalia",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_mammalia)

png(file = paste0("results/plots/O1_0_110ky_", core_target, ".png"), width = 20, height = 13, units = "in", res = 600)
egg::ggarrange(
  plotBacillariophyta, plotChlorophyta, plotHaptophyta, plotOthers,
  plotArthropoda, plotActinopteri, plotMammalia, nrow = 1)
dev.off()










# Set parameters for the plots
x_max_s <- 150000
x_min_s <- 0
core_target  <- c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")
trans_applied <- "sqrt"
min_n <- 5
min_reads <- 1000


# Filter table based on samples and taxa
samples <- species_stats_mapping %>%
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    core %in% core_target
  ) %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))

taxa <- species_stats_mapping %>%
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    core %in% core_target
  ) %>%
  group_by(species) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))



datafilt <- species_stats_mapping %>%
  mutate(n_reads = ifelse(species == "Emiliania huxleyi", n_reads/100, n_reads)) %>%
  # Filter data
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+25000,
    species %in% taxa$species[taxa$n_reads >= min_reads & taxa$n >= min_n],
    # label %in% samples$label[samples$n_reads >= 500 & samples$n >= 3],
    core %in% core_target
  ) %>% 
  # Complete the data
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  # Add taxa and core information
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  left_join(cdata)
 

## Primary producers
plotBacillariophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Bacillariophyta",
  trans = trans_applied, y_labels = TRUE, 
  core_levels = core_target, title = "Bacillariophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Bacillariophyta)

plotChlorophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Chlorophyta",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Chlorophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Chlorophyta)

plotHaptophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Haptophyta",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Haptophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Haptophyta)

plotOthers <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = c("Pelagophyceae", "Bolidophyceae", "Cryptophyceae", "Dictyochophyceae", "Eustigmatophyceae", "Polycystinea"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Other primary producers",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_opp)

## Intermediary and top consumers
plotArthropoda <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = c("Arthropoda", "Foraminifera", "Cercozoa", "Priapulida"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Microorg. consummers",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_arthopoda)

plotActinopteri <- plot_reads_by_taxa(datafilt,
  taxa_col = "class", taxa_value = c("Actinopteri", "Chondrichthyes"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Fishes",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_fishes)

plotMammalia <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Mammalia",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Mammalia",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_mammalia)


png(file = paste0("results/plots/O1_0_110ky.png"), width = 45, height = 13, units = "in", res = 600)
egg::ggarrange(
  plotBacillariophyta, plotChlorophyta, plotHaptophyta, plotOthers,
  plotArthropoda, plotActinopteri, plotMammalia, nrow = 1)
dev.off()




















# Set parameters for the plots
x_max_s <- 500000
x_min_s <- 0
core_target  <- c("GeoB25202_R1")
trans_applied <- "sqrt"
min_n <- 1
min_reads <- 1000


# Filter table based on samples and taxa
samples <- species_stats_mapping %>%
  filter(
    fit == "good" | y_bp < 10000,
    # fit == "good",
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    core %in% core_target
  ) %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))

taxa <- species_stats_mapping %>%
  filter(
    fit == "good" | y_bp < 10000,
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    core %in% core_target
  ) %>%
  group_by(species) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))

datafilt <- species_stats_mapping %>%
  # Filter data
  filter(
    fit == "good" | y_bp < 10000,
    y_bp > x_min_s-25000 & y_bp < x_max_s+25000,
    species %in% taxa$species[taxa$n_reads >= min_reads & taxa$n >= min_n],
    # label %in% samples$label[samples$n_reads >= 500 & samples$n >= 3],
    core %in% core_target
  ) %>%
  # Complete the data
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  # Add taxa and core information
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  left_join(cdata)

## Primary producers
plotBacillariophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Bacillariophyta",
  trans = trans_applied, y_labels = TRUE, 
  core_levels = core_target, title = "Bacillariophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Bacillariophyta)

plotChlorophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Chlorophyta",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Chlorophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Chlorophyta)

plotHaptophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Haptophyta",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Haptophyta",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_Haptophyta)

plotOthers <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = c("Pelagophyceae", "Bolidophyceae", "Cryptophyceae", "Dictyochophyceae", "Eustigmatophyceae", "Polycystinea"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Other primary producers",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_opp)

## Intermediary and top consumers
plotArthropoda <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = c("Arthropoda", "Foraminifera", "Cercozoa", "Priapulida"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Microorg. consummers",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_arthopoda)

plotActinopteri <- plot_reads_by_taxa(datafilt,
  taxa_col = "class", taxa_value = c("Actinopteri", "Chondrichthyes"),
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Fishes",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_fishes)

plotMammalia <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Mammalia",
  trans = trans_applied, y_labels = FALSE, 
  core_levels = core_target, title = "Mammalia",
  x_min = x_min_s, x_max = x_max_s, ncols = 1, palette = colors_mammalia)

png(file = paste0("results/plots/O3_0_500ky_", core_target, ".png"), width = 20, height = 13, units = "in", res = 600)
egg::ggarrange(
  plotBacillariophyta, 
  plotChlorophyta, 
  plotHaptophyta, 
  plotOthers,
  plotArthropoda,
  plotActinopteri, 
  plotMammalia, nrow = 1)
dev.off()



species_stats_mapping %>%
  filter(species == "Boreogadus saida") %>%
  arrange(desc(n_reads)) %>%
  View()

genus_stats_mapping %>%
  filter(genus == "Boreogadus") %>%
  arrange(desc(n_reads)) %>%
  ggplot(., aes(x = y_bp, y = n_reads)) +
  geom_point()+
  facet_nested(.~core, scales = "free")+
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
  coord_flip()











detail <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = c("Pelagophyceae", "Polycystinea", "Dinophyceae", "Bolidophyceae", "Priapulimorpha", "Priapulimorpha", "Echinoidea", "Thecofilosea"),
  trans = NULL, y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s, facet = TRUE)

png(file = "test_detail.png", width = 40, height = 20, units = "in", res = 600)
detail + theme(aspect.ratio = 4/1, legend.position = "bottom")
dev.off()

"
Look at the serie covering from 0 to ~120 ky bp:
  




Look at the serie covergin from ~90ky bp to ~350ky bp, potentially capturing three interglaciar stages and two glaciar stages:

Some taxa was presetn all across the record:
  Fish: Gradus morhua (Atlantic cod)
  Coscinodiscophyceae: Thalassiosira rotula.
  Chlorophyta: Bathycoccus prasinos. 

At the starting of the observed peaks in terms of total number of reads and observed species, there is a recurren appearence of
Clupea harengus (Atlantic herring). Together with Clupea harengus, members of the genus Sebastes (Scorpaenidae) are also present, 
but with a lower robustnes. While Clupea harengus appears, the number of reads of top consumers (Mammalia) remains low and only
Balaenoptera acutorostrata (Common minke whale) and, in lower numbers, Balaenoptera musculus (Blue whale) are present. Regarding 
Bacillariophyceae, Farilariopsis cylindrus, the most abundant member of the group, reach its minimum abundances. In some events,
Pseudo-nitzschia delicatissima and Pseudo-nitzschia cuspidata are also present together with Clupea harengus, but reach maximum 
abundances at later stages, when Clupea harengus is not present. Thalassiosira rotula remains as the most abundant 
Coscinodiscophyceae, and Micromonas comoda appears at some HP events, together with the other Chlorophyta species Bathycoccus 
prasinos. 






Fishes:
  Gadus morhua (Atlantic cod)
  Clupea harengus (Atlantic herring)
  Sebastes umbrosus (Rougheye rockfish)
  Boreogadus saida (Arctic cod)
  Hippoglossus hippoglossus (Atlantic halibut)
  Hyppoglossus stenolepis (Pacific halibut)
  Tarulus bubalis (Atlantic wolffish)

Mammals:
  Balaenoptera acutorostrata (Common minke whale)
  Balaenoptera musculus (Blue whale)
  Halichoerus grypus (Grey seal)
  Phoca vitulina (Harbour seal)
  Orcinus orca (Orca)
  Physeter catodon (Sperm whale)
  Globicephala melas (Long-finned pilot whale)
  Mirounga angustirostris (Northern elephant seal)
  Leptonychotes weddellii (Weddell seal)
  Neophocaena schauinslandi (Hawaiian monk seal)
  Monodon monoceros (Narwhal)
  Mirounga leonina (Southern elephant seal)


Gadur morhua is present acrros the whole record including both glacial and interglacial periods.

Taxa appreareing together with Clupea harengus:



"

plotMammalia <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Mammalia",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s, facet = TRUE)














# Filtrado de especies válidas
samples <- genus_stats_mapping %>%
  filter(fit == "good") %>%
  group_by(label) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))

taxa <- genus_stats_mapping %>%
  filter(fit == "good") %>%
  group_by(genus) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  arrange(desc(n))

# View range
x_max_s <- 120000
x_min_s <- 0

datafilt <- genus_stats_mapping %>%
  # Filter data
  filter(
    y_bp > x_min_s-25000 & y_bp < x_max_s+20000,
    genus %in% taxa$genus[taxa$n_reads >= 1000 & taxa$n >= 10],
    # label %in% samples$label[samples$n_reads >= 500 & samples$n >= 5],
    core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")
  ) %>%
  # Complete the data
  group_by(label, genus) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        genus = unique(genus),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  # Add taxa and core information
  left_join(genus_stats_mapping %>% select(phylum, class, order, family, genus) %>% distinct()) %>%
  left_join(cdata)

nreads <- datafilt %>%
  group_by(label, core, y_bp) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  mutate(core = fct_relevel(core, c("ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
  ggplot(aes(x = y_bp, y = n_reads)) +
			geom_vline(
			data = datafilt %>% 
				filter(grepl("Clupea", genus) & n_reads > 0) %>%
				mutate(core = fct_relevel(core, c("ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
				select(label, core, y_bp),
				aes(xintercept = y_bp), linetype = 1, color = "grey", linewidth = 2, alpha = 0.5) +
			geom_point()+
      geom_line() +
			facet_nested(~core, scales = "free") +
			guides(fill = guide_legend(ncol = 1)) +
			scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
			coord_flip(xlim = c(x_max_s, x_min_s))+
      # coord_flip()+
			theme(aspect.ratio = 4/1, legend.position = "bottom")

nsp <- datafilt %>%
  group_by(label, core, y_bp) %>%
  filter(n_reads > 0) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  mutate(core = fct_relevel(core, c("ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
  ggplot(aes(x = y_bp, y = n)) +
			geom_vline(
			data = datafilt %>% 
				filter(grepl("Clupea", genus) & n_reads > 0) %>%
				mutate(core = fct_relevel(core, c("ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
				select(label, core, y_bp),
				aes(xintercept = y_bp), linetype = 1, color = "grey", linewidth = 2, alpha = 0.5) +
			geom_point()+
      geom_line() +
			facet_nested(~core, scales = "free") +
			guides(fill = guide_legend(ncol = 1)) +
			scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
			coord_flip(xlim = c(x_max_s, x_min_s))+
      # coord_flip()+
			theme(aspect.ratio = 4/1, legend.position = "bottom")



## Intermediary and top consumers
plotArthropoda <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Arthropoda",
  trans = "sqrt", y_labels = TRUE, 
  core_levels = c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotActinopteri <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Actinopteri",
  trans = sqrt, y_labels = FALSE, 
  core_levels = c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotMammalia <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Mammalia",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")


png(file = "test.png", width = 20, height = 10, units = "in", res = 600)
egg::ggarrange(plotArthropoda, plotActinopteri, plotMammalia, nrow = 1)
dev.off()



plotActinopteri <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Actinopteri",
  trans = "sqrt", y_labels = TRUE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotMammalia <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Mammalia",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotBacillariophyceae <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Bacillariophyceae",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotCoscinodiscophyceae <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = "Coscinodiscophyceae",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotChlorophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Chlorophyta",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotHaptophyta <- plot_reads_by_taxa(datafilt, 
  taxa_col = "phylum", taxa_value = "Haptophyta",
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

plotOther <- plot_reads_by_taxa(datafilt, 
  taxa_col = "class", taxa_value = c("Pelagophyceae", "Polycystinea", "Dinophyceae", "Bolidophyceae", "Priapulimorpha", "Priapulimorpha", "Echinoidea", "Thecofilosea"),
  trans = "sqrt", y_labels = FALSE, 
  core_levels = c("ST13", "GeoB25202_R1", "GeoB25202_R2"), 
  x_min = x_min_s, x_max = x_max_s,
  resolution = "genus")

png(file = "test.png", width = 40, height = 10, units = "in", res = 600)
egg::ggarrange(nsp, nreads, plotActinopteri, plotMammalia, plotBacillariophyceae, plotCoscinodiscophyceae, 
  plotChlorophyta, plotHaptophyta, plotOther, nrow = 1)
dev.off()


































plotA <- species_stats_mapping %>%
  filter(species %in% taxa$species[taxa$n_reads >= 1000 & taxa$n >= 5]) %>%
  # filter(y_bp > 200000 & y_bp < 380000) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(class == "Mammalia") %>%
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(species) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  left_join(cdata) %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(species = fct_reorder(species, -sum_n_reads)) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = species)) +
    geom_vline(
      data = species_stats_mapping %>% 
        filter(species == "Clupea harengus", core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>% mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>% select(label, core, y_bp),
      aes(xintercept = y_bp), linetype = 1, color = "grey", linewidth = 2, alpha = 0.3)+
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus) +

    # scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4), limits = c(380000, 180000))+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4), limits = c(380000, 190000))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")+
    theme(aspect.ratio = 3/1)
plotA





## Reads by group
plotA <- species_stats_mapping %>%
  mutate(group = case_when(
    species %in% fishes ~ "Fishes",
    species %in% mammals ~ "Mammals",
    species %in% fitoplankton ~ "Phytoplankton",
    species %in% zooplankton ~ "Zooplankton",
    TRUE ~ "Other"
  )) %>% #select(phylum, class, order, family, genus, species, group) %>% distinct() %>% View()
  group_by(label, group) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  left_join(cdata) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(aes(x = y_bp, y = n_reads)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core*group, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")
plotA




## Reads by group
plotA <- species_stats_mapping %>%
  mutate(group = case_when(
    species %in% fishes ~ "Fishes",
    species %in% mammals ~ "Mammals",
    species %in% fitoplankton ~ "Phytoplankton",
    species %in% zooplankton ~ "Zooplankton",
    TRUE ~ "Other"
  )) %>% #select(phylum, class, order, family, genus, species, group) %>% distinct() %>% View()
  filter(group == "Phytoplankton") %>%
  group_by(species, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata) %>%
  group_by(species) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(species = fct_reorder(species, -sum_n_reads)) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = species)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_manual(values = paired_genus) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")
plotA

## Reads by group
plotA <- species_stats_mapping %>%
  filter(fit == "good") %>%
  mutate(group = case_when(
    species %in% fishes ~ "Fishes",
    species %in% mammals ~ "Mammals",
    species %in% fitoplankton ~ "Phytoplankton",
    species %in% zooplankton ~ "Zooplankton",
    TRUE ~ "Other"
  )) %>% #select(phylum, class, order, family, genus, species, group) %>% distinct() %>% View()
  filter(group == "Phytoplankton") %>%
  group_by(label) %>%
  summarise(n_reads = n(), .groups = "drop") %>%
  left_join(cdata) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(aes(x = y_bp, y = n_reads)) +
    geom_point()+
    geom_line() +
    # geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    # geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_manual(values = paired_genus) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")
plotA


## Reads by group
plotA <- species_stats_mapping %>%
  filter(core == "GeoB25202_R1") %>%
  filter(fit == "good") %>%
  mutate(group = case_when(
    species %in% fishes ~ "Fishes",
    species %in% mammals ~ "Mammals",
    species %in% fitoplankton ~ "Phytoplankton",
    species %in% zooplankton ~ "Zooplankton",
    TRUE ~ "Other"
  )) %>% #select(phylum, class, order, family, genus, species, group) %>% distinct() %>% View()
  # filter(group == "Phytoplankton") %>%
  group_by(label, group) %>%
  summarise(n_reads = n(), .groups = "drop") %>%
  complete(.,
        group = unique(group),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(aes(x = y_bp, y = n_reads, color = group)) +
    geom_point()+
    geom_line() +
    # geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    # geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Total number of reads assigned at species level", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_manual(values = paired_genus) +
    facet_nested(~core*group, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    # scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")
plotA











# Librerías necesarias
library(corrplot)

target_cores <- c("ST13", "GeoB25202_R1", "GeoB25202_R2")
# target_core <- "ST13"
results <- list()

for (target_core in target_cores) {

  data <- species_stats_mapping %>%
    filter(species %in% fitoplankton) %>%
    filter(y_bp > 245000 & y_bp < 350000) %>%
    filter(grepl(target_core, core)) %>%
    filter(fit == "good")

  # Filtrado de muestras válidas
  samples <- data %>%
    group_by(label) %>%
    summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
    arrange(n_reads)

  # Filtrado de especies válidas
  taxa <- data %>%
    group_by(species) %>%
    summarise(n_reads = mean(n_reads), n = n(), .groups = "drop") %>%
    arrange(desc(n))

  # Umbrales mínimos
  min_reads_sample <- 1
  min_sp_sample <- 1
  min_reads_specie <- 1
  min_sp_specie <- 1

  # Datos en formato ancho
  data_wide <- species_stats_mapping %>%
    group_by(label) %>%
    mutate(n_reads = sum(n_reads)) %>%
    ungroup() %>%
    filter(label %in% as.character(samples %>% filter(n_reads >= min_reads_sample & n >= min_sp_sample) %>% pull(label))) %>%
    filter(species %in% as.character(taxa %>% filter(n_reads >= min_reads_specie & n >= min_sp_specie) %>% pull(species))) %>%
    filter(fit == "good") %>%
    select(label, species, n_reads) %>%
    pivot_wider(names_from = species, values_from = n_reads, values_fill = 0) %>%
    column_to_rownames("label")

  # Función para calcular matriz de p-valores
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    return(p.mat)
  }

  # Matrices de correlación y p-valores
  M <- cor(data_wide, method = "spearman", use = "pairwise.complete.obs")
  p.mat <- cor.mtest(data_wide, method = "spearman")

  # Guardar resultados
  results[[paste0("M_", target_core)]] <- M
  results[[paste0("p.mat_", target_core)]] <- p.mat
}


# # Visualización: solo correlaciones significativas (p < 0.05)
png(file = "test.png", width = 15, height = 7, units = "in", res = 600)
par(mfrow = c(1, 3))
# ST8 corrplot
# corrplot(results$M_ST8, method = "color", order = "hclust", type = "upper", 
#          p.mat = results$p.mat_ST8, sig.level = 0.05, insig = "blank", main = "ST8")
# ST13 corrplot
corrplot(results$M_ST13, method = "color", order = "hclust", type = "upper", 
         p.mat = results$p.mat_ST13, sig.level = 0.05, insig = "blank", main = "ST13")
# GeoB25202_R1 corrplot
corrplot(results$M_GeoB25202_R1, method = "color", order = "hclust", type = "upper", 
         p.mat = results$p.mat_GeoB25202_R1, sig.level = 0.05, insig = "blank", main = "GeoB25202 (R1)")
# # GeoB25202_R2 corrplot
corrplot(results$M_GeoB25202_R2, method = "color", order = "hclust", type = "upper", 
         p.mat = results$p.mat_GeoB25202_R2, sig.level = 0.05, insig = "blank", main = "GeoB25202 (R2)")
# Resetea layout a uno solo
par(mfrow = c(1, 1))
dev.off()

corrplot(results$M_GeoB25202_R2, method = "color", order = "hclust", type = "upper", 
         p.mat = results$p.mat_GeoB25202_R2, sig.level = 0.05, insig = "blank", main = "GeoB25202 (R1)")

# corrplot(M, 
#   # method = 'square', 
#   type = "upper", 
#   diag = TRUE, 
#   order = 'hclust',
#   p.mat = p.mat, 
#   sig.level = 0.05,
#   # addrect = 5,
#   insig = "blank", 
#   rect.col = 'blue',
#   rect.lwd = 3)
















data <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1")) %>%
  filter(fit == "good") %>%
  filter(y_bp > 50000) %>%
  filter(y_bp < 110000)

abundance <- data %>%
  group_by(species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>% 
  arrange(desc(n_reads))

abundance %>% print(n = 70)

plotA <- data %>%
  filter(species %in% c(c(abundance$species[1:18], "Fragilariopsis cylindrus"))) %>%
  # filter(family == "Noelaerhabdaceae" | species %in% c("Triparma laevis", "Sebastes umbrosus", "Pseudo-nitzschia cuspidata", "Thalassiosira rotula", "Balaenoptera musculus", "Gadus morhua", "Phoca vitulina", "Clupea harengus", "Balaenoptera acutorostrata", "Halichoerus grypus")) %>%
  # mutate(taxa = ifelse(family == "Noelaerhabdaceae", family, species)) %>% 
  mutate(taxa = species) %>% 
  # mutate(taxa = ifelse(is.na(taxa), "Other", taxa)) %>%
  group_by(taxa, label, fit) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>% #View()
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata %>% select(label, core, y_bp)) %>%
  group_by(taxa) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  # mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  # mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  # mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(taxa = fct_reorder(taxa, -sum_n_reads)) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "identity", alpha = 1, width = 1000) +
    geom_area(aes(color = taxa), stat = "identity", position = "identity", alpha = 0.3) +
    # scale_fill_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "brown2"), labels = c("Noelaerhabdaceae" = "Family Noelaerhabdaceae", "Other" = "Other taxa")) +
    # scale_color_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "brown2"), guide = "none") +
    labs(title = "Reads assigned to family Noelaerhabdaceae", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(core~taxa, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    # scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")

png(file = "test.png", width = 20, height = 18, units = "in", res = 600)
plotA
dev.off()


















## Noelaerhabdaceae
plotA <- species_stats_mapping %>%
  mutate(taxa = ifelse(family == "Noelaerhabdaceae", "Noelaerhabdaceae", "Other")) %>% 
  mutate(taxa = ifelse(is.na(taxa), "Other", taxa)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata %>% select(label, core, y_bp)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(taxa = fct_relevel(taxa, "Other", "Noelaerhabdaceae")) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2) +
    geom_area(aes(color = taxa), stat = "identity", position = "stack", alpha = 0.5) +
    scale_fill_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "darkgrey"), labels = c("Noelaerhabdaceae" = "Family Noelaerhabdaceae", "Other" = "Other taxa")) +
    scale_color_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "darkgrey"), guide = "none") +
    labs(title = "Reads assigned to family Noelaerhabdaceae", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    # scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")+
    theme(aspect.ratio = 3/1)



# png(file = "results/eukaryal/plots/prelim/02_dominance_noelaerhabdaceae.png", width = 12, height = 8, units = "in", res = 600)
plotA
# dev.off()



## Herring
plotB <- species_stats_mapping %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  filter(taxa %in% c("Clupea harengus")) %>% 
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1000)+
    geom_area(aes(color = taxa), stat = "identity", position = "stack", alpha = 0.3) +
    labs(title = "Reads assigned to Clupea harengus (Atlantic herring)", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = "brown2")+
    scale_color_manual(values = "brown2", guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)
    
plotB

# png(file = "results/eukaryal/plots/prelim/04_herring.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plotA + theme(legend.position = "none"), plotB + theme(legend.position = "none"), nrow = 1, widths = c(1, 1))
# dev.off()



abundance <- species_stats_mapping %>%
  filter(core %in% c("ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good") %>%
  filter(y_bp < 250000) %>%
  filter(y_bp > 180000) %>%
  group_by(species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>% 
  arrange(desc(n_reads)) #%>% 
  # print(n = 30)

species_stats_mapping %>%
  filter(fit == "good") %>%
  group_by(phylum, class, order, family, genus, species) %>%
  summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") %>%
  filter(n > 1) %>%
  # select(-n_reads, -n) %>%
  arrange() %>% 
  print(n = 170)

plotA <- species_stats_mapping %>%
  filter(core %in% c("ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good") %>%
  filter(y_bp < 250000) %>%
  filter(y_bp > 180000) %>%
  filter(species %in% c(abundance$species[1:18])) %>%
  # filter(family == "Noelaerhabdaceae" | species %in% c("Triparma laevis", "Sebastes umbrosus", "Pseudo-nitzschia cuspidata", "Thalassiosira rotula", "Balaenoptera musculus", "Gadus morhua", "Phoca vitulina", "Clupea harengus", "Balaenoptera acutorostrata", "Halichoerus grypus")) %>%
  # mutate(taxa = ifelse(family == "Noelaerhabdaceae", family, species)) %>% 
  mutate(taxa = species) %>% 
  # mutate(taxa = ifelse(is.na(taxa), "Other", taxa)) %>%
  group_by(taxa, label, fit) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>% #View()
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata %>% select(label, core, y_bp)) %>%
  group_by(taxa) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  # mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  # mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  # mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(taxa = fct_reorder(taxa, -sum_n_reads)) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "identity", alpha = 1, width = 1000) +
    geom_area(aes(color = taxa), stat = "identity", position = "identity", alpha = 0.3) +
    # scale_fill_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "brown2"), labels = c("Noelaerhabdaceae" = "Family Noelaerhabdaceae", "Other" = "Other taxa")) +
    # scale_color_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "brown2"), guide = "none") +
    labs(title = "Reads assigned to family Noelaerhabdaceae", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(core~taxa, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
    scale_y_sqrt(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    # scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")

plotA





abundance <- genus_stats_mapping %>%
  filter(core == "GeoB25202_R2") %>%
  filter(fit == "good") %>%
  filter(y_bp < 400000) %>%
  filter(y_bp > 250000) %>%
  group_by(genus) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>% 
  arrange(desc(n_reads)) %>% 
  filter(genus != "Emiliania")
  # print(n = 30)

plotA <- genus_stats_mapping %>%
  filter(core == "GeoB25202_R2") %>%
  # filter(fit == "good") %>%
  filter(y_bp < 400000) %>%
  filter(y_bp > 250000) %>%
  filter(genus %in% c(abundance$genus[1:18])) %>%
  # filter(family == "Noelaerhabdaceae" | species %in% c("Triparma laevis", "Sebastes umbrosus", "Pseudo-nitzschia cuspidata", "Thalassiosira rotula", "Balaenoptera musculus", "Gadus morhua", "Phoca vitulina", "Clupea harengus", "Balaenoptera acutorostrata", "Halichoerus grypus")) %>%
  # mutate(taxa = ifelse(family == "Noelaerhabdaceae", family, species)) %>% 
  mutate(taxa = genus) %>% 
  # mutate(taxa = ifelse(is.na(taxa), "Other", taxa)) %>%
  group_by(taxa, label, fit) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>% #View()
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata %>% select(label, core, y_bp)) %>%
  group_by(taxa) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  # mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  # mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  # mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(taxa = fct_reorder(taxa, -sum_n_reads)) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "identity", alpha = 1, width = 1000) +
    geom_area(aes(color = taxa), stat = "identity", position = "identity", alpha = 0.3) +
    # scale_fill_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "brown2"), labels = c("Noelaerhabdaceae" = "Family Noelaerhabdaceae", "Other" = "Other taxa")) +
    # scale_color_manual(values = c("Noelaerhabdaceae" = "darkgreen", "Other" = "brown2"), guide = "none") +
    labs(title = "Reads assigned to family Noelaerhabdaceae", x = "", y = "Number of reads") +
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core*taxa, scales = "free") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    # scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip() +
    labs(fill = "")
plotA










## Fragilaropsis
plotA <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  filter(taxa %in% c("Fragilariopsis cylindrus")) %>% 
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads assigned to Fragilaropsis cylindrus", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = "darkcyan")+
    scale_color_manual(values = "darkcyan", guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "results/eukaryal/plots/prelim/05_fragilariopsis_cylindrus.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 5))
dev.off()



## Balaenoptera
plotA <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(species %in% c("Balaenoptera acutorostrata", "Balaenoptera musculus")) %>% 
  rename(taxa = species) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads assigned to Balaenoptera", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus)+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "results/eukaryal/plots/prelim/03_balaenoptera.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 5))
dev.off()


## Calanus
plotA <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  filter(taxa %in% c("Calanus glacialis", "Calanus hyperboreus", "Calanus finmarchicus")) %>% 
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "identity", alpha = 0.5) +
    labs(title = "Reads assigned to Calanus finmarchicus, glacialis and hyperboreus", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus)+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "results/eukaryal/plots/prelim/06_calanus.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 5))
dev.off()



## Calanus 50 reads
species_stats_mapping_50 <- read.csv("results/eukaryal/data_processing/species50_stats_more1k_dmg_and_non-dmg.csv") %>%
	filter(!grepl("GeoB25206", core)) %>%
  left_join(cdata %>% select(label, mis, temp))


plotA <- species_stats_mapping_50 %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  filter(taxa %in% c("Calanus glacialis", "Calanus hyperboreus", "Calanus finmarchicus", "Metridia gerlachei")) %>% 
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "identity", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "identity", alpha = 0.5) +
    labs(title = "Reads assigned to Calanus finmarchicus, glacialis and hyperboreus (50 reads)", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus)+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "results/eukaryal/plots/prelim/06_calanus_down-to-50.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 5))
dev.off()



## Phocidae
plotA <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  filter(taxa %in% c("Phoca vitulina", "Halichoerus grypus", "Mirounga angustirostris", "Leptonychotes weddellii", "Neomonachus schauinslandi", "Mirounga leonina")) %>% 
  group_by(taxa) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(taxa = fct_reorder(taxa, -sum_n_reads)) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads assigned to family Phocidae", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus)+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "results/eukaryal/plots/prelim/07_phocidae.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 5))
dev.off()



## Clupea and Fragilaropsis
plotA <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  filter(taxa %in% c("Fragilariopsis cylindrus", "Clupea harengus")) %>% 
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "identity", alpha = 0.5) +
    labs(title = "Reads assigned to Fragilaropsis cylindrus and Clupea harengus", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = c("brown2", "darkcyan"))+
    scale_color_manual(values = c("brown2", "darkcyan"), guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    scale_y_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

png(file = "results/eukaryal/plots/prelim/06_fragilariopsis_and_clupea.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 10))
dev.off()



## Miscelaneous
plotA <- species_stats_mapping %>%
  filter(core %in% c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  rename(taxa = species) %>%
  select(taxa, label, n_reads) %>%
  tidyr::complete(.,
        taxa = unique(taxa),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  group_by(taxa, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  # filter(taxa %in% c("Gadus morhua")) %>% 
  # filter(taxa %in% c("Thalassiosira rotula")) %>% 
  # filter(taxa %in% c("Phoca vitulina", "Halichoerus grypus")) %>% 
  # filter(taxa %in% c("Phoca vitulina", "Halichoerus grypus", "Mirounga angustirostris", "Leptonychotes weddellii", "Neomonachus schauinslandi", "Mirounga leonina")) %>% 
  filter(taxa %in% c("Bathycoccus prasinos", "Thalassiosira rotula", "Gadus morhua", "Phaeocystis antarctica", "Triparma laevis")) %>% 
  group_by(taxa) %>%
  mutate(sum_n_reads = sum(n_reads)) %>%
  ungroup() %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  mutate(taxa = fct_reorder(taxa, sum_n_reads)) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(aes(color = taxa), stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads assigned to Calanus finmarchicus, glacialis and hyperboreus", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core*taxa, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus)+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
    # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# png(file = "results/eukaryal/plots/prelim/06_calanus.png", width = 12, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, nrow = 1, widths = c(1, 1, 10))
# dev.off()














library(microbiome)
library(phyloseq)
library(robCompositions)
library(vegan)
library(gghighlight)
library(ggh4x)
library(ggrepel)
library(indicspecies)


##### Species level #####
# Load species data and modify it 
sp_ps <- readRDS(file = "results/eukaryal/taxonomy/ps-filt-sp.rds") %>%
    prune_samples(!grepl("ST5", sample_data(.)$core), .)

sample_occurrence <- sp_ps %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    filter(Abundance > 0) %>%
    group_by(Sample, y_bp, core) %>%
    summarise(n = n(), sum_reads = sum(Abundance), .groups = "drop") %>%
    arrange(desc(n))

sp_occurrence <- sp_ps %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    filter(Abundance > 0) %>%
    group_by(species) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(n)

species_ps_filt <- sp_ps %>%
    prune_samples(sample_names(.) %in% as.character(sample_occurrence %>% filter(n >= 5 & sum_reads > 1000) %>% pull(Sample)), .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_taxa(taxa_names(.) %in% as.character(sp_occurrence %>% filter(n >= 10) %>% pull(species)), .)

species_ps_filt <- prune_taxa(taxa_names(species_ps_filt) != "Emiliania huxleyi", species_ps_filt)
sp_trans <- microbiome::transform(species_ps_filt, "clr")

# PCoA
ps <- sp_trans

# Original data
data <- species_ps_filt %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    rename(taxa = OTU, n_reads = Abundance)

# Distance matrix
dist <- ps %>%
    otu_table() %>%
    t() %>%
    data.frame() %>%
    vegan::vegdist(., method = "euclidean")

#### PCoA ####
pcoa_result <- cmdscale(dist, eig = TRUE, k = 2)

# Obtener los valores propios (eigenvalues)
eigenvalues <- pcoa_result$eig

# Calcular el porcentaje de varianza explicada por cada eje
varianza_explicada <- eigenvalues / sum(eigenvalues)

# Crear un data frame con las coordenadas y los labels
coordenadas <- as.data.frame(pcoa_result$points) %>%
    rename_at(vars(c("V1", "V2")), ~ c("PC1", "PC2")) %>% 
    rownames_to_column("label") %>%
    inner_join(cdata, by = "label") %>%
    arrange(label)

# Taxa contribution
data_c <- ps %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    rename(taxa = OTU, n_reads = Abundance) %>%
    filter(n_reads != 0) %>%
    select(taxa, label, n_reads) %>%
    arrange(label) %>%
    pivot_wider(names_from = label, values_from = n_reads, values_fill = 0) %>%
    column_to_rownames("taxa") %>%
    as.data.frame() %>%
    t() 
# rownames(data_c) == coordenadas$label

# Calcular loadings como correlaciones entre las abundancias de las especies y las coordenadas PCoA
loadings <- cor(data_c, pcoa_result$points)

# Crear un data frame para los loadings
loadings_df <- as.data.frame(loadings)
loadings_df$variable <- rownames(loadings_df)


# Calcular loadings como correlaciones entre las variables originales y las coordenadas PCoA
contribuciones <- apply(abs(loadings), 1, sum)
contribuciones_ordenadas <- sort(contribuciones, decreasing = TRUE)
contribuciones_ordenadas %>%
    as.data.frame() %>%
    rownames_to_column("taxa") %>%
    rename("contribution" = 2) %>%
    mutate(taxa = fct_reorder(taxa, contribution)) %>%
    ggplot(., aes(x = taxa, y = contribution))+
        geom_point()+
        coord_flip()

# Crear un data frame para los loadings
loadings_df <- as.data.frame(loadings) %>%
    rownames_to_column("variable") %>%
    rename_at(vars(c("V1", "V2")), ~ c("PC1", "PC2")) %>%
    filter(complete.cases(PC1) & complete.cases(PC2)) %>%
    mutate(contribution = abs(PC1) + abs(PC2)) %>%
    arrange(desc(contribution)) %>%
    left_join(data %>% select(species, genus, family, order, class, phylum) %>% distinct(), by = c("variable" = "species")) #%>%
    # mutate(PC1 = PC1*10,
    #        PC2 = PC2*10)

# Añadir flechas al gráfico
# Clustering
hc <- hclust(dist(as.data.frame(pcoa_result$points)), method = "ward.D2")
dist <- dist(as.data.frame(pcoa_result$points))

k_values <- 2:10
avg_sil_values <- sapply(k_values, function(k) avg_silhouette(hc, k, dist))

optimal_k <- k_values[which.max(avg_sil_values)]
optimal_k
# optimal_k <- 3
plot(k_values, avg_sil_values, xlab="k",ylab="av. silhouette",type="b", pch=19)

clusters_opt <- cutree(hc, k = optimal_k) %>%
    as.data.frame() %>% 
    rownames_to_column("label") %>% 
    rename_at(vars(2), ~ c("cluster_opt")) 

# Plot the dendrogram to decide number of clusters
dendro_data <- ggdendro::dendro_data(hc)

dendro_metadata <- data.frame(pcoa_result$points) %>%
    rownames_to_column("label") %>%
    left_join(cdata) %>%
    select(label, core, y_bp)

label_data <- merge(dendro_data$labels, dendro_metadata)
clusters <- cutree(hc, k = optimal_k)
label_data$cluster <- clusters[as.character(label_data$label)]

segment_data <- dendro_data$segments
segment_data$cluster <- label_data$cluster[match(segment_data$xend, label_data$x)]

dendroPlot <- ggplot() +
    geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(data = label_data, aes(x = x, y = y-5, color = as.factor(cluster)), size = 1) +
    # scale_color_manual(values = c("2"="#F8766D", "1"="#00BFC4", "3"="#09037f", "4"="#0d9d19"))+
    # scale_color_manual(values = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00"))+
    labs(x = "", y = "Height", color = "Core")+
    # ylim(c(-10, 10))+
    theme_minimal()+
    theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grob_secundario <- ggplotGrob(dendroPlot)

# Join cluster information to nmds points
coordenadas_clust <- coordenadas %>%
    left_join(clusters_opt, by = "label")
    
coordPlot <- ggplot(coordenadas_clust, aes(x = PC1, y = PC2)) +

    geom_point(aes(fill = as.factor(cluster_opt), size = derep_reads, shape = core))+ 
    # scale_fill_manual(values = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00"))+

    stat_ellipse(data = coordenadas_clust %>% filter(cluster_opt == 1), level = 0.95, type = "norm", linetype = 2)+
    stat_ellipse(data = coordenadas_clust %>% filter(cluster_opt == 2), level = 0.95, type = "norm", linetype = 2)+
    stat_ellipse(data = coordenadas_clust %>% filter(cluster_opt == 3), level = 0.95, type = "norm", linetype = 2)+
    stat_ellipse(data = coordenadas_clust %>% filter(cluster_opt == 4), level = 0.95, type = "norm", linetype = 2)+

    # geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2, color = phylum), linewidth = 1.5, alpha = 0.2)+
    # geom_text(data = loadings_df %>% filter(contribution > quantile(contribution, .8)) %>% mutate(phylum = ifelse(is.na(phylum), "Other", phylum)), aes(label = variable, color = phylum), check_overlap = TRUE)+

    scale_size_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    scale_shape_manual(values = c(21, 22, 24, 25, 26))+
    guides(fill = "none")+
    labs(x = paste("PCoA1 (", round(varianza_explicada[1] * 100, 2), "% explained variance)", sep = ""),
        y = paste("PCoA2 (", round(varianza_explicada[2] * 100, 2), "% explained variance)", sep = ""),
        title = "PCoA clustering of samples based on (clr-transformed) species composition",
        fill = "Cluster",
        shape = "Core",
        color = "Phylum",
        size = "Number of\nunique reads")

coordPlotMod <- coordPlot +
  annotation_custom(grob_secundario, xmin = -10, xmax = 0, ymin = 10.5, ymax = 5.5)+
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1)), 
    xmin = -10, xmax = -0.1, ymin = 10.5, ymax = 6)+
    theme(aspect.ratio = 1/1)


png(file = "results/eukaryal/plots/prelim/07_pcoa.png", width = 10, height = 10, units = "in", res = 600)
coordPlotMod
dev.off()


#  Indicator species analysis
data_ind <- ps %>%
    otu_table() %>%
    as.data.frame() %>%
    t()  %>%
    as.data.frame() %>%
    rownames_to_column(var = "label") %>% 
    full_join(clusters_opt) %>%
    column_to_rownames("label")

data_ind_obs <- data_ind %>% select(-cluster_opt)
data_ind_obs_pos <- data_ind_obs+ data_ind %>% select(-cluster_opt) %>% min()*-1

indval <- multipatt(data_ind_obs_pos, data_ind$cluster_opt, func = "IndVal.g", duleg = TRUE)

# summary(indval, alpha = 0.01)

indv_group <- indval$sign %>%
    as.data.frame() %>%
    filter(p.value <= 0.01) %>%
    rownames_to_column(var = "taxa") %>%
    select(taxa, index) %>%
    mutate(index = paste0("Group_", index))

core_ranges <- cdata %>% group_by(core) %>% summarise(min_y_bp = min(y_bp), max_y_bp = max(y_bp)) %>% ungroup()

species_reads <- sp_ps %>%
    otu_table() %>%
    as.data.frame() %>%
    rownames_to_column(var = "taxa") %>%
    as_tibble() %>%
    pivot_longer(cols = -taxa, names_to = "label", values_to = "n_reads") %>%
    left_join(sp_ps %>% sample_data() %>% as.data.frame() %>% as_tibble()) %>%
    left_join(sp_ps %>% tax_table() %>% as.data.frame() %>% rownames_to_column(var = "taxa")) %>%
    group_by(species) %>%
    mutate(total_reads = sum(n_reads),
        occurrence = sum(n_reads > 0)) %>%
    ungroup() %>%
    mutate(relative_reads = total_reads/sum(total_reads)) %>%
    mutate(species = fct_reorder(species, total_reads)) %>%
    mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    # mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    left_join(cdata %>% group_by(core) %>% summarise(min_y_bp = min(y_bp), max_y_bp = max(y_bp)) %>% ungroup()) %>%
    group_by(taxa) %>%
    left_join(indv_group) %>%
    filter(!is.na(index)) %>%
    mutate(species = fct_reorder(species, total_reads)) %>%
    left_join(clusters_opt) %>%
    # mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))
    mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(core = str_replace(core, "25202_R1", " (R1)")) %>%
    mutate(core = str_replace(core, "25202_R2", " (R2)")) %>%
    mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB (R1)", "GeoB (R2)"))

# Set color palette
tax_paths_unique <- species_reads %>% pull(species) %>% unique()
colors_genus <- paired_genus[1:length(tax_paths_unique)]
names(colors_genus) <- tax_paths_unique


samples <- species_reads %>%
    filter(index == paste0("Group_", cluster_opt)) %>%
    mutate(xmin = y_bp - sqrt(y_bp)*3,
           xmax = y_bp + sqrt(y_bp)*3,
           ymin = 0,
           ymax = max(n_reads, na.rm = TRUE)) %>%
    select(label, core, index, cluster_opt, xmin, xmax, ymin, ymax) %>%
    distinct()

# ggplot(species_reads, aes(x = y_bp, y = n_reads, fill = species))+
#     geom_rect(aes(xmin = min(min_y_bp), xmax = min_y_bp, ymin = 0, ymax = Inf), fill = "#e4e2e2", alpha = 0.1) + 
#     geom_rect(aes(xmin = max_y_bp, xmax = max(max_y_bp), ymin = 0, ymax = Inf), fill = "#e4e2e2", alpha = 0.1) +
#     geom_rect(data = samples, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf), 
#         fill = "#1F78B4", alpha = 0.2, inherit.aes = FALSE) +
#     geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1)+
#     geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
#     facet_nested(index~., scales = "free")+
#     guides(color = guide_legend(override.aes = list(shape = c(16, 17, 18, 19))))+
#     scale_fill_manual(values = colors_genus, guide = guide_legend(reverse = TRUE, ncol = 1))+
#     scale_y_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
#     scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000))+
#     labs(x = "Age", y = "Number of reads", fill = "", title = "Group 1 (clustering without Emiliania)")


# Filtrar los datos para los puntos que quieres resaltar
samples_1 <- species_reads %>% 
    ungroup() %>%
    filter(cluster_opt == 1, index == "Group_1") %>%
    mutate(xmin = y_bp - sqrt(y_bp)*3,
           xmax = y_bp + sqrt(y_bp)*3,
           ymin = 0,
           ymax = max(n_reads, na.rm = TRUE)) %>%
    select(label, core, index, cluster_opt, xmin, xmax, ymin, ymax) %>%
    distinct()

plotA <- ggplot(species_reads %>% filter(index == "Group_1"), aes(x = y_bp, y = n_reads, fill = species))+
            geom_rect(data = samples_1, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf), 
              fill = "#1F78B4", alpha = 0.2, inherit.aes = FALSE) +
            geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1)+
            geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
            facet_nested(.~core, scales = "free")+
            scale_fill_manual(values = colors_genus, guide = guide_legend(reverse = TRUE, ncol = 1))+
            scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 2))+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
            # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
            labs(x = "", y = "", fill = "", title = "Samples clustering and indicator species abundance suggesting three different regimes")+
            coord_flip()+
            theme(aspect.ratio = 4/1, legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())


samples_2 <- species_reads %>% 
    ungroup() %>%
    filter(cluster_opt == 2, index == "Group_2") %>%
    mutate(xmin = y_bp - sqrt(y_bp)*3,
           xmax = y_bp + sqrt(y_bp)*3,
           ymin = 0,
           ymax = max(n_reads, na.rm = TRUE)) %>%
    select(label, core, index, cluster_opt, xmin, xmax, ymin, ymax) %>%
    distinct()

plotB <- ggplot(species_reads %>% filter(index == "Group_2"), aes(x = y_bp, y = n_reads, fill = species))+
            geom_rect(data = samples_2, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf), 
              fill = "#33A02C", alpha = 0.2, inherit.aes = FALSE) +
            geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1)+
            geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
            facet_nested(.~core, scales = "free")+
            scale_fill_manual(values = colors_genus, guide = guide_legend(reverse = TRUE, ncol = 1))+
            scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 2))+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
            # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
            labs(x = "", y = "Number of reads", fill = "", title = "")+
            coord_flip()+
            theme(aspect.ratio = 4/1, legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())


samples_3 <- species_reads %>% 
    ungroup() %>%
    filter(cluster_opt == 3, index == "Group_3") %>%
    mutate(xmin = y_bp - sqrt(y_bp)*3,
           xmax = y_bp + sqrt(y_bp)*3,
           ymin = 0,
           ymax = max(n_reads, na.rm = TRUE)) %>%
    select(label, core, index, cluster_opt, xmin, xmax, ymin, ymax) %>%
    distinct()

plotC <- ggplot(species_reads %>% filter(index == "Group_3"), aes(x = y_bp, y = n_reads, fill = species))+
            geom_rect(data = samples_3, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf), 
              fill = "#E31A1C", alpha = 0.2, inherit.aes = FALSE) +
            geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1)+
            geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
            facet_nested(.~core, scales = "free")+
            scale_fill_manual(values = colors_genus, guide = guide_legend(reverse = TRUE, ncol = 1))+
            scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 2))+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
            # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
            labs(x = "", y = "", fill = "", title = "")+
            coord_flip()+
            theme(aspect.ratio = 6/1, legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())


png(file = "results/eukaryal/plots/prelim/08_clusters.png", width = 13, height = 8, units = "in", res = 600)
egg::ggarrange(plottemp, plotmi, plotA, plotB, plotC, nrow = 1, widths = c(1,1,3,3,3))
dev.off()
































png(file = "results/eukaryal/plots/prelim/number_reads_good_fit.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  group_by(label, fit) |>
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  group_by(label) %>%
  mutate(relative_n_reads = n_reads / sum(n_reads)) %>%
  ungroup() %>%
  select(label, fit, n_reads) %>%
  tidyr::complete(.,
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  inner_join(cdata) |> 
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = fit)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Number of reads classified as damaged and non-damaged at a species level", x = "Age (y bp)", y = "Number of reads") +
    guides(fill = guide_legend(ncol = 1)) +
    facet_nested(~core, scales = "free")+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")
dev.off()

species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, fit, n_reads) %>%
  tidyr::complete(.,
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  group_by(label, fit) |>
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  inner_join(cdata) |> 
  group_by(core, fit) %>%
  summarise(mean_reads = round(mean(n_reads),0), median_reads = round(median(n_reads), 0), .groups = "drop") %>%
  arrange(fit, core) %>%
  knitr::kable()
  


png(file = "results/eukaryal/plots/prelim/damage_good_fit.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  ggplot(aes(x = y_bp, y = A_b, color = fit, size = n_reads)) +
    geom_point()+
    labs(title = "Damage estimates for species classified as damaged and non-damaged at a species level", x = "Age (y bp)", y = "Damage (A_b)", color = "", size = "Number of\nreads") +
    guides(fill = guide_legend(ncol = 1)) +
    facet_nested(~core, scales = "fixed")+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    # scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")
dev.off()


png(file = "results/eukaryal/plots/prelim/number_species_good_fit.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  group_by(label, fit) |>
  count() |>
  ungroup() |>
  select(label, fit, n) %>%
  tidyr::complete(.,
        label = unique(label),
        fit = unique(fit),
        fill = list(n = 0)) %>%
  inner_join(cdata) |> 
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  ggplot(aes(x = y_bp, n, fill = fit)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Number of species classified as damaged and non-damaged", x = "Age (y bp)", y = "Number of species") +
    guides(fill = guide_legend(ncol = 1)) +
    facet_nested(~core, scales = "free")+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    coord_flip()+
    labs(fill = "")
dev.off()



nsp_fit <- species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  group_by(label) %>%
  count() %>%
  ungroup()

png(file = "results/eukaryal/plots/prelim/weighted_number_species_good_fit.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  group_by(label, fit) |>
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  group_by(label) %>%
  mutate(relative_n_reads = n_reads / sum(n_reads)) %>%
  ungroup() %>% 
  left_join(nsp_fit) %>%
  mutate(proportion = n*relative_n_reads) %>%
  # group_by(label) %>%
  # mutate(sum_proportion = sum(proportion)) %>%
  # ungroup() %>%
  # View()
  select(label, fit, proportion) %>%
  tidyr::complete(.,
        label = unique(label),
        fit = unique(fit),
        fill = list(proportion = 0)) %>%
  inner_join(cdata) |> 
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  ggplot(aes(x = y_bp, y = proportion, fill = fit)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads weighted number of species classified as damaged and non-damaged", x = "Age (y bp)", y = "Number of species") +
    guides(fill = guide_legend(ncol = 1)) +
    facet_nested(~core, scales = "free")+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    coord_flip()+
    labs(fill = "")
dev.off()

species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  group_by(label, fit) %>%
  count() %>%
  ungroup() %>%
  left_join(cdata) %>%
  group_by(core, fit) %>%
  summarise(mean_sp = round(mean(n),0), .groups = "drop") %>%
  arrange(fit, core) %>%
  knitr::kable()








sample_stats <- species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  group_by(label) %>%
  count() %>%
  ungroup() %>%
  full_join(readswoE) %>% 
  left_join(cdata %>% select(label, core, y_bp))
  
ps_unfilt <- readRDS(file = "results/eukaryal/taxonomy/phyloseq_unfiltered-sp.rds") %>% 
  prune_samples(!grepl("GeoB25206", sample_data(.)$core), .) #%>% 
  # prune_samples(sample_names(.) %in% as.character(sample_stats %>% filter(n_reads > 1000 & n > 5) %>% pull(label)), .)


ps <- ps_unfilt %>%
  # prune_taxa(taxa_names(.) != "Emiliania huxleyi", .) %>%
  # microbiome::transform("compositional")
  microbiome::transform("clr")

ps %>%
  speedyseq::psmelt() %>%
  as_tibble() %>%
  rename(taxa = OTU, n_reads = Abundance) %>%
  filter(species == targets[12]) %>%
  select(species, n_reads, label, core, y_bp) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(aes(x = y_bp, y = n_reads))+
    geom_point()+
    geom_line()+
    facet_nested(~core, scales = "free")+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    coord_flip()
  








targets <- species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  left_join(cdata %>% select(label, core)) %>%
  select(core, species, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        core = unique(core),
        fill = list(n_reads = 0)) %>% 
  group_by(core, species) %>%
  summarise(mean_reads = mean(n_reads), sum_reads = sum(n_reads), .groups = "drop") %>%
  group_by(species) %>%
  mutate(total_reads = sum(sum_reads)) %>%
  ungroup() %>% 
  filter(total_reads > quantile(total_reads, 0.9))  %>%
  pull(species) %>%
  unique()



target <- targets[12]
target <- "Emiliania huxleyi"
print(target)  
plot <- species_stats_mapping %>% 
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  # filter(fit == "Damaged") |>
  select(label, species, fit, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  filter(species == target) %>% 
  left_join(cdata) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = fit)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = paste0("Reads assigned to ", target), x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")

png(file = paste0("results/eukaryal/plots/prelim/main_taxa/", janitor::make_clean_names(target), ".png"), width = 15, height = 8, units = "in", res = 600)
plot
dev.off()

target <- targets[12]
target <- c("Phaeocystis antarctica", "Clupea harengus", "Fragilariopsis cylindrus", "Chaetoceros gelidus")
plot <- species_stats_mapping %>% 
  filter(!grepl("GeoB25206", core)) %>%
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  filter(species %in% target) %>% 
  left_join(cdata) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  filter(fit == "Damaged") |>
  ggplot(aes(x = y_bp, y = n_reads, fill = species)) +
    geom_area(aes(color = species), stat = "identity", position = "stack", alpha = 0.7) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    # geom_point(data = . %>% filter(fit != "Damaged"), aes(color = fit), alpha = 1)+
    labs(title = paste0("Reads assigned to ", target), x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = c("Phaeocystis antarctica" = "#1f77b4", 
                             "Clupea harengus" = "darkred", 
                             "Fragilariopsis cylindrus" = "#2ca02c", 
                             "Chaetoceros gelidus" = "#9933CC"))+
    scale_color_manual(values = c("Phaeocystis antarctica" = "#1f77b4", 
                             "Clupea harengus" = "darkred", 
                             "Fragilariopsis cylindrus" = "#2ca02c", 
                             "Chaetoceros gelidus" = "#9933CC"))+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")

plot





# plot <- species_stats_mapping %>% 
#   filter(!grepl("GeoB25206", core)) |>
#   mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
#   # filter(fit == "Damaged") |>
#   select(label, species, fit, n_reads) %>%
#   tidyr::complete(.,
#         species = unique(species),
#         label = unique(label),
#         fit = unique(fit),
#         fill = list(n_reads = 1)) %>%
#   filter(species == target) %>% 
#   left_join(cdata %>% select(label, core, y_bp)) |>
#   group_by(core) %>%
#   arrange(-y_bp) %>%
#   mutate(log2fold = log2(n_reads/lag(n_reads))) %>%
#   ungroup() %>%
#   mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#   mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
#   filter(fit == "Damaged") %>%
#   # filter(log2fold > 0) %>%
#   ggplot(aes(x = y_bp, y = log2fold, fill = fit)) +
#     geom_point()+
#     geom_line()+
#     # geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
#     # geom_area(stat = "identity", position = "stack", alpha = 0.5) +
#     labs(title = paste0("Reads assigned to ", target), x = "Age (y bp)", y = "Number of reads")+
#     guides(fill = guide_legend(nrow = 1)) +
#     facet_nested(~core, scales = "free")+
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
#     scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#     coord_flip()+
#     labs(fill = "")
# plot



mi <- read.table(file = "/projects/fernandezguerra/people/ngm902/ROCS/associated_data/Lisiecki2005_copy.txt", sep = "\t", header = T)
plotmi <- mi %>%
  mutate(y_bp = Time_ka * 1e3) %>%
  filter(y_bp <= 497443) %>%
  ggplot(., aes(x = y_bp, y = Benthic_d18O))+
    geom_line()+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_reverse()+
    coord_flip()+
    labs(x = "Age (y bp)", y = "Benthic d18O", title = "Marine Isotope Stage")

png(file = "results/eukaryal/plots/prelim/mis.png", width = 3, height = 8, units = "in", res = 600)
plotmi
dev.off()



species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  filter(core != "ST5") %>%
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  filter(species %in% targets[15:20]) %>% 
  left_join(cdata %>% select(label, core, y_bp)) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  filter(n_reads > 0) %>%
  # head()
  ggplot(aes(x = y_bp, y = n_reads, shape = fit, color = core)) +
    geom_point()+
    labs(title = paste0("Reads assigned to ", target), x = "Age (y bp)", y = "Number of reads")+
    # geom_smooth(method = "loess", se = FALSE, span = 0.3)+
    facet_nested(~species, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")















png(file = "results/eukaryal/plots/prelim/emiliania_species_dmg.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  filter(grepl("Emiliania huxleyi", species)) %>% 
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  left_join(cdata) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = fit)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads assigned to Emiliania huxleyi", x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")
dev.off()

png(file = "results/eukaryal/plots/prelim/emiliania_species.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, n_reads) %>%
  # filter(grepl("Emiliania huxleyi", species)) %>% 
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(cdata) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(species = fct_relevel(species, "Other", "Emiliania huxleyi")) %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = species)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.7) +
    scale_fill_manual(values = c("Emiliania huxleyi" = "orangered2", "Other" = "grey"))+
    labs(title = "Reads assigned to Emiliania huxleyi", x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")
dev.off()


species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, n_reads) %>%
  # filter(grepl("Emiliania huxleyi", species)) %>% 
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania_huxleyi", "Other")) %>%
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  pivot_wider(names_from = species, values_from = n_reads) %>%
  mutate(all = Emiliania_huxleyi + Other) %>%
  left_join(cdata %>% select(label, core, y_bp)) %>%
  ggplot(., aes(x = all, y = Emiliania_huxleyi))+
    geom_abline(xintercept = 0, slope = 1, linetype = "dashed", color = "red", alpha = .5)+
    geom_point()+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    facet_nested(~core, scales = "fixed")+
    theme(aspect.ratio = 1)


species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  # filter(fit == "Damaged") |> 
  select(label, species, n_reads) %>% 
  # filter(grepl("Emiliania huxleyi", species)) %>% 
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  group_by(label, species) %>% 
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  # tidyr::complete(.,
  #       species = unique(species),
  #       label = unique(label),
  #       fill = list(n_reads = 0)) %>%
  pivot_wider(names_from = species, values_from = n_reads) %>% 
  rename(emiliania = 2) %>% 
  mutate(emiliania = ifelse(is.na(emiliania), 0, emiliania)) %>% 
  mutate(Other = ifelse(is.na(Other), 0, Other)) %>% 
  mutate(total = emiliania+Other, proportion = emiliania/total) %>% 
  left_join(cdata %>% select(label, core)) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = core, y = proportion))+
    geom_lv(alpha = 0.5)+
    geom_jitter(width = 0.1, height = 0)+
    coord_flip()+
    labs(y = "Proportion of reads assigned to Emiliania huxleyi", x = "")



png(file = "results/eukaryal/plots/prelim/non-emiliania_species.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  group_by(label, species, fit) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  left_join(cdata) |>
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  filter(fit == "Damaged") %>%
  # head()
  ggplot(aes(x = y_bp, y = n_reads, fill = fit)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    labs(title = "Reads assigned to other than Emiliania huxleyi", x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "fixed")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(10000, 100000, 400000))+
    coord_flip()+
    labs(fill = "")
dev.off()

readswoE <- species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  group_by(label, species, fit) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  filter(fit == "Damaged") %>%
  select(-fit, -species) 

png(file = "results/eukaryal/plots/prelim/species_vs_reads_woEmiliania.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  group_by(label) %>%
  count() %>%
  ungroup() %>%
  full_join(readswoE) %>% 
  left_join(cdata %>% select(label, core, y_bp)) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = n_reads, y = n))+
    geom_point()+
    facet_nested(~core, scales = "fixed")+
    scale_x_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1e3, 1e4, 5e4, 1e5, 2e5, 3e5))+
    theme(aspect.ratio = 1)+
    labs(x = "Number of reads assigned to other than Emiliania huxleyi", y = "Number of species classified as damaged")
dev.off()

species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  group_by(label) %>%
  count() %>%
  ungroup() %>%
  full_join(readswoE) %>% 
  left_join(cdata %>% select(label, core, y_bp)) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(pass = ifelse(n_reads > 1e4 & n > 10, "yes", "no")) %>%
  group_by(core) %>% 
  summarise(total = n(), more_de_10k = sum(pass == "yes"), .groups = "drop") %>%
  mutate(proportion = less_de_10k/total) %>%
  head()
  head()
















species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  group_by(label, species, fit) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  left_join(cdata %>% select(label, core, y_bp)) %>%
  filter(fit == "Damaged") %>%
  ggplot(aes(x = n_reads, y = core)) +
      geom_violin()+
      geom_point()+
      facet_nested(core~., scales = "free_y")+
      scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
      



species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(core, label, species, n_reads) %>%
  mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>%
  # filter(fit == "Damaged") %>%
  group_by(label, species,  core) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  group_by(core) %>%
  summarise(median_reads = round(median(n_reads),0), mean_reads = round(mean(n_reads),0), .groups = "drop") %>%
  knitr::kable()


png(file = "results/eukaryal/plots/prelim/occurrence_species.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  # mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>% 
  left_join(cdata %>% select(label, core)) %>%
  group_by(species, core) %>%
  summarise(
    mean_reads = mean(n_reads),
    mean_reads_wp = round(mean(n_reads[n_reads > 0]), 0),
    occurrences = sum(n_reads > 0),
    sum_reads = sum(n_reads),
    .groups = "drop") %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = mean_reads_wp, y = occurrences))+
    ggrepel::geom_text_repel(data = . %>% filter(occurrences > 50 | mean_reads > 10e2), aes(label = species, color = species), nudge_x = 50, nudge_y = 1, size = 3)+
    # scale_color_manual(values = paired_genus)+
    geom_point()+
    facet_nested(~core, scales = "fixed")+
    theme(aspect.ratio = 1, legend.position = "none")+
    scale_x_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    labs(x = "Mean number of reads per sample", y = "Number of samples where present")
dev.off()



species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>% 
  left_join(cdata %>% select(label, core)) %>%
  group_by(species, core) %>%
  summarise(
    mean_reads = mean(n_reads),
    mean_reads_wp = round(mean(n_reads[n_reads > 0]), 0),
    occurrences = sum(n_reads > 0),
    sum_reads = sum(n_reads),
    .groups = "drop") %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%  
  head()


species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>% 
  left_join(cdata %>% select(label, core)) %>%
  group_by(species, core) %>%
  summarise(
    mean_reads = mean(n_reads),
    mean_reads_wp = round(mean(n_reads[n_reads > 0]), 0),
    occurrences = sum(n_reads > 0),
    sum_reads = sum(n_reads),
    .groups = "drop") %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, genus, species) %>% distinct()) %>%
  filter_all(any_vars(grepl("Chaetoceros", .))) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  arrange(desc(mean_reads)) %>%
  print(n = 20)
  ggplot(., aes(x = mean_reads, y = occurrences))+
    ggrepel::geom_text_repel(data = . %>% filter(occurrences > 20 | mean_reads > 10e2), aes(label = species, color = species), nudge_x = 50, nudge_y = 1, size = 3)+
    geom_point()+
    facet_nested(~core, scales = "fixed")+
    theme(aspect.ratio = 1, legend.position = "none")+
    scale_x_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    labs(x = "Mean number of reads per sample", y = "Number of samples where present")

target <- "Calanus glacialis|Calanus finmarchicus|Calanus hyperboreus"
plotA <- species_stats_mapping %>% 
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, genus, species) %>% distinct()) %>%
  filter_all(any_vars(grepl(target, .))) %>%
  left_join(cdata) |>
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  mutate(species = fct_reorder(species, -sum_reads)) %>%
  filter(fit == "Damaged") %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = species)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    scale_fill_manual(values = paired_genus)+
    labs(title = paste0("Reads assigned to ", target), x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    # scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")

# png(file = paste0("results/eukaryal/plots/prelim/main_taxa/", janitor::make_clean_names(target), ".png"), width = 15, height = 8, units = "in", res = 600)
png(file = paste0("results/eukaryal/plots/prelim/main_taxa/", "calanus", ".png"), width = 15, height = 8, units = "in", res = 600)
plotA
dev.off()



target <- "Balaenoptera acutorostrata|Balaenoptera musculus|Orcinus orca"
plotA <- species_stats_mapping %>% 
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  select(label, species, fit, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fit = unique(fit),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, genus, species) %>% distinct()) %>%
  filter_all(any_vars(grepl(target, .))) %>%
  left_join(cdata) |>
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(fit = fct_relevel(fit, "Non-damaged", "Damaged")) %>%
  mutate(species = fct_reorder(species, -sum_reads)) %>%
  filter(fit == "Damaged") %>%
  ggplot(aes(x = y_bp, y = n_reads, fill = species)) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 2)+
    geom_area(stat = "identity", position = "stack", alpha = 0.5) +
    scale_fill_manual(values = paired_genus)+
    labs(title = paste0("Reads assigned to ", target), x = "Age (y bp)", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    # scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
    scale_y_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    coord_flip()+
    labs(fill = "")

png(file = paste0("results/eukaryal/plots/prelim/main_taxa/", "balaenoptera_orcinus", ".png"), width = 15, height = 8, units = "in", res = 600)
plotA
dev.off()




species_stats_mapping %>% filter_all(any_vars(grepl("Orcinus", .))) %>% View()
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>% 
  left_join(cdata %>% select(label, core)) %>%
  group_by(species, core) %>%
  summarise(
    mean_reads = mean(n_reads),
    mean_reads_wp = round(mean(n_reads[n_reads > 0]), 0),
    occurrences = sum(n_reads > 0),
    sum_reads = sum(n_reads),
    .groups = "drop") %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, genus, species) %>% distinct()) %>%
  filter_all(any_vars(grepl("Cetacea", .))) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  arrange(desc(mean_reads)) %>%
  print(n = 20)
  ggplot(., aes(x = mean_reads, y = occurrences))+
    ggrepel::geom_text_repel(data = . %>% filter(occurrences > 20 | mean_reads > 10e2), aes(label = species, color = species), nudge_x = 50, nudge_y = 1, size = 3)+
    geom_point()+
    facet_nested(~core, scales = "fixed")+
    theme(aspect.ratio = 1, legend.position = "none")+
    scale_x_continuous(trans = "sqrt", labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    labs(x = "Mean number of reads per sample", y = "Number of samples where present")











species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  # mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>% 
  left_join(cdata %>% select(label, core)) %>%
  group_by(species, core) %>%
  summarise(
    mean_reads = mean(n_reads),
    mean_reads_wp = round(mean(n_reads[n_reads > 0]), 0),
    occurrences = sum(n_reads > 0),
    sum_reads = sum(n_reads),
    .groups = "drop") %>%
  select(species, core, sum_reads) %>%
  distinct() %>%
  group_by(core) %>%
  mutate(relative_reads = sum_reads/sum(sum_reads)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -relative_reads)) %>%
  ggplot(., aes(x = core, y = relative_reads, fill = species))+
    geom_bar(stat = "identity")+
    # coord_flip()+
    # scale_fill_manual(values = paired_genus)+
    labs(y = "Relative number of reads", x = "Core")+
    theme(aspect.ratio = 1, legend.position = "none")
dev.off()





png(file = "results/eukaryal/plots/prelim/top10percent_species.png", width = 15, height = 8, units = "in", res = 600)
species_stats_mapping %>%
  filter(!grepl("GeoB25206", core)) |>
  mutate(fit = ifelse(fit == "good", "Damaged", "Non-damaged")) |>
  filter(fit == "Damaged") %>%
  select(label, species, n_reads) %>%
  # mutate(species = ifelse(species == "Emiliania huxleyi", "Emiliania huxleyi", "Other")) %>%
  filter(!grepl("Emiliania huxleyi", species)) %>% 
  group_by(label, species) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  left_join(cdata %>% select(label, core)) %>%
  # arrange(desc(sum_reads)) %>%
  select(core, species, n_reads) %>%
  tidyr::complete(.,
        species = unique(species),
        core = unique(core),
        fill = list(n_reads = 0)) %>% 
  group_by(core, species) %>%
  summarise(mean_reads = mean(n_reads), sum_reads = sum(n_reads), .groups = "drop") %>%
  group_by(species) %>%
  mutate(total_reads = sum(sum_reads)) %>%
  ungroup() %>% 
  filter(total_reads > quantile(total_reads, 0.9)) %>%
  mutate(species = fct_reorder(species, total_reads)) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = species, y = mean_reads))+
    geom_bar(stat = "identity")+
    facet_nested(~core)+
    coord_flip()+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    labs(y = "Mean sample reads", x = "Top species by number of reads", title = "Mean number of reads for top 10% species by number of reads. Only damaged reads are considered.")
dev.off()


# agg_species_stats_mapping_dmg %>%
#   write.csv(., "/projects/fernandezguerra/people/ngm902/ROCS/results/euk/data_processing/species_stats_raw.csv", row.names = FALSE)

# # Filter based on mapping, damage and reference length
# species_mapping_dmg_filtered <- agg_species_stats_mapping_dmg %>%
#   filter(rm == "keep", fit == "good", sum_reference_length > 1000)

# species_mapping_dmg_filtered %>%
#   write.csv(., "/projects/fernandezguerra/people/ngm902/ROCS/results/euk/data_processing/species_stats_more1k.csv", row.names = FALSE)

# species_mapping_all_filtered <- agg_species_stats_mapping_dmg %>%
#   filter(rm == "keep", sum_reference_length > 1000)

# species_mapping_all_filtered %>%
#   write.csv(., "/projects/fernandezguerra/people/ngm902/ROCS/results/euk/data_processing/species_stats_more1k_dmg_and_non-dmg.csv", row.names = FALSE)

  

















species_stats_mapping %>%
  select(phylum, class, order) %>%
  distinct() %>%
  filter(grepl("Actinopteri", class)) %>%
  # filter(grepl("Lamniformes", order)) %>%
  head()
  

fishes <- c("Gadiformes", "Salmoniformes", "Clupeiformes", "Perciformes", "Osmeriformes",  "Aulopiformes", "Rajiformes", "Cyprinodontiformes", "Cichliformes", "Spariformes", "Centrarchiformes", "Carangiformes",
  "Lamniformes",        "Anguilliformes")


core_ranges <- cdata %>% filter(!grepl("GeoB25206", core)) %>% group_by(core) %>% summarise(min_y_bp = min(y_bp), max_y_bp = max(y_bp)) %>% ungroup()

plotA <- species_stats_mapping %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(y_bp < 30000) %>%
  group_by(species, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  filter(order %in% c("Lamniformes")) %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  left_join(core_ranges) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = species)) +
    geom_rect(aes(xmin = 0, xmax = min_y_bp, ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) + 
    geom_rect(aes(xmin = max_y_bp, xmax = max(max_y_bp), ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) +
    geom_area(aes(color = species), stat = "identity", position = "stack", alpha = .5) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 200)+
    labs(title = "Reads assigned to Carcharodon carcharias", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus, guide = "none")+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(limits = c(30000, 0), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 30000, 1000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")+
    theme(aspect.ratio = 4/1, legend.position = "none")

png(file = "results/eukaryal/plots/prelim/fishes/great_white_shark.png", width = 10, height = 8, units = "in", res = 600)
plotA
dev.off()


plotA <- species_stats_mapping %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(y_bp < 30000) %>%
  group_by(species, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  filter(class == "Actinopteri") %>%
  filter(order == "Gadiformes") %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  left_join(core_ranges) %>%
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -sum_reads)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = species)) +
    geom_rect(aes(xmin = 0, xmax = min_y_bp, ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) + 
    geom_rect(aes(xmin = max_y_bp, xmax = max(max_y_bp), ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) +
    geom_area(aes(color = species), stat = "identity", position = "stack", alpha = .5) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 200)+
    labs(title = "Reads assigned to order Gadiformes", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus, guide = "none")+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(limits = c(30000, 0), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 30000, 1000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")
    # theme(aspect.ratio = 4/1, legend.position = "none")

png(file = "results/eukaryal/plots/prelim/fishes/Gadiformes.png", width = 10, height = 8, units = "in", res = 600)
plotA
dev.off()


plotA <- species_stats_mapping %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(y_bp < 30000) %>%
  group_by(species, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  filter(class == "Actinopteri") %>%
  filter(order == "Clupeiformes") %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  left_join(core_ranges) %>%
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -sum_reads)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = species)) +
    geom_rect(aes(xmin = 0, xmax = min_y_bp, ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) + 
    geom_rect(aes(xmin = max_y_bp, xmax = max(max_y_bp), ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) +
    geom_area(aes(color = species), stat = "identity", position = "stack", alpha = .5) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 200)+
    labs(title = "Reads assigned to order Clupeiformes", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus, guide = "none")+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(limits = c(30000, 0), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 30000, 1000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")

png(file = "results/eukaryal/plots/prelim/fishes/Clupeiformes.png", width = 10, height = 8, units = "in", res = 600)
plotA
dev.off()


plotA <- species_stats_mapping %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(y_bp < 30000) %>%
  group_by(species, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  filter(class == "Actinopteri") %>%
  filter(order == "Perciformes") %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  left_join(core_ranges) %>%
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -sum_reads)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = species)) +
    geom_rect(aes(xmin = 0, xmax = min_y_bp, ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) + 
    geom_rect(aes(xmin = max_y_bp, xmax = max(max_y_bp), ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) +
    geom_area(aes(color = species), stat = "identity", position = "stack", alpha = .5) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 200)+
    labs(title = "Reads assigned to order Perciformes", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus, guide = "none")+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(limits = c(30000, 0), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 30000, 1000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")

png(file = "results/eukaryal/plots/prelim/fishes/Perciformes.png", width = 10, height = 8, units = "in", res = 600)
plotA
dev.off()


plotA <- species_stats_mapping %>%
  filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  filter(fit == "good" | y_bp < 10000) %>%
  filter(y_bp < 30000) %>%
  group_by(species, label) %>%
  summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  tidyr::complete(.,
        species = unique(species),
        label = unique(label),
        fill = list(n_reads = 0)) %>%
  left_join(species_stats_mapping %>% select(phylum, class, order, family, genus, species) %>% distinct()) %>%
  filter(class == "Actinopteri") %>%
  filter(order == "Osmeriformes") %>%
  left_join(cdata %>% select(label, core, y_bp, mis)) %>%
  left_join(core_ranges) %>%
  group_by(species) %>%
  mutate(sum_reads = sum(n_reads)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, -sum_reads)) %>%
  mutate(core = str_replace(core, "_R1", " (R1)")) %>%
  mutate(core = str_replace(core, "_R2", " (R2)")) %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202 (R1)", "GeoB25202 (R2)")) %>%
  ggplot(., aes(x = y_bp, y = n_reads, fill = species)) +
    geom_rect(aes(xmin = 0, xmax = min_y_bp, ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) + 
    geom_rect(aes(xmin = max_y_bp, xmax = max(max_y_bp), ymin = 0, ymax = Inf), fill = "#D1D3D5", alpha = 0.3) +
    geom_area(aes(color = species), stat = "identity", position = "stack", alpha = .5) +
    geom_bar(stat = "identity", position = "stack", alpha = 1, width = 200)+
    labs(title = "Reads assigned to order Osmeriformes", x = "", y = "Number of reads")+
    guides(fill = guide_legend(nrow = 1)) +
    facet_nested(~core, scales = "free")+
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = paired_genus, guide = "none")+
    scale_color_manual(values = paired_genus, guide = "none")+
    scale_x_reverse(limits = c(30000, 0), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 30000, 1000))+
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::pretty_breaks(n = 3))+
    coord_flip()+
    labs(fill = "")

png(file = "results/eukaryal/plots/prelim/fishes/Osmeriformes.png", width = 10, height = 8, units = "in", res = 600)
plotA
dev.off()


