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
taxonomy <- tax_data %>% select(domain, phylum, class, order, family, genus, species) %>% distinct()
# reference_summary <- read_tsv("results/microbial/damage/damage-classification-depositional/reference-summary.tsv")
# sample_baselines <- read_tsv("results/microbial/damage/damage-classification-depositional/sample-baselines.tsv")


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





# ============================
# Damaged-only descriptive exploration (preliminary)
# ============================

dmg_only_plot_dir <- "plots/microbial/preliminary"
dmg_only_res_dir <- "results/microbial/damage/damage-classification-depositional/exploration"
dir.create(dmg_only_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dmg_only_res_dir, recursive = TRUE, showWarnings = FALSE)

dmg_only_data <- tax_data %>%
  filter(is_dmg == "Damaged") %>%
  mutate(
    domain = ifelse(is.na(domain) | domain == "", "Unclassified_domain", domain),
    phylum = ifelse(is.na(phylum) | phylum == "", "Unclassified_phylum", phylum),
    class = ifelse(is.na(class) | class == "", "Unclassified_class", class),
    order = ifelse(is.na(order) | order == "", "Unclassified_order", order),
    family = ifelse(is.na(family) | family == "", "Unclassified_family", family),
    genus = ifelse(is.na(genus) | genus == "", "Unclassified_genus", genus),
    species = ifelse(is.na(species) | species == "", "Unclassified_species", species)
  )

dmg_only_meta <- cdata %>%
  select(label, core, y_bp) %>%
  distinct()

dmg_only_stats_path <- file.path(dmg_only_res_dir, "sample-classification-stats.tsv")
dmg_only_stats <- if (file.exists(dmg_only_stats_path)) {
  read_tsv(dmg_only_stats_path, show_col_types = FALSE) %>%
    select(label, pct_reads_final_damaged)
} else {
  tibble(label = character(), pct_reads_final_damaged = numeric())
}

dmg_only_sample_summary <- dmg_only_data %>%
  group_by(label) %>%
  summarise(
    total_damaged_abundance = sum(abundance, na.rm = TRUE),
    n_species_damaged = n_distinct(species[!is.na(species) & species != "" & !grepl("^Unclassified", species)]),
    n_genus_damaged = n_distinct(genus[!is.na(genus) & genus != "" & !grepl("^Unclassified", genus)]),
    n_family_damaged = n_distinct(family[!is.na(family) & family != "" & !grepl("^Unclassified", family)]),
    .groups = "drop"
  ) %>%
  left_join(dmg_only_meta, by = "label") %>%
  left_join(dmg_only_stats, by = "label") %>%
  select(label, core, y_bp, total_damaged_abundance, n_species_damaged, n_genus_damaged, n_family_damaged, pct_reads_final_damaged) %>%
  arrange(core, y_bp)

write_tsv(dmg_only_sample_summary, file.path(dmg_only_res_dir, "damaged_only_sample_summary.tsv"))

dmg_only_n_samples <- n_distinct(dmg_only_sample_summary$label)

dmg_only_top_taxa <- dmg_only_data %>%
  group_by(domain, phylum, class, order, family, genus, species) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    n_samples_present = n_distinct(label[abundance > 0]),
    .groups = "drop"
  ) %>%
  mutate(
    prevalence_fraction = n_samples_present / ifelse(dmg_only_n_samples > 0, dmg_only_n_samples, NA_real_)
  ) %>%
  arrange(desc(total_abundance), desc(prevalence_fraction), desc(n_samples_present))

write_tsv(dmg_only_top_taxa, file.path(dmg_only_res_dir, "damaged_only_top_taxa.tsv"))

dmg_only_dom_phylum <- dmg_only_data %>%
  group_by(core, phylum) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(core) %>%
  mutate(rank = dense_rank(desc(abundance))) %>%
  filter(rank <= 3) %>%
  arrange(core, rank) %>%
  summarise(dominant_phyla = paste0(phylum, collapse = "; "), .groups = "drop")

dmg_only_dom_genus <- dmg_only_data %>%
  group_by(core, genus) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(core) %>%
  mutate(rank = dense_rank(desc(abundance))) %>%
  filter(rank <= 3) %>%
  arrange(core, rank) %>%
  summarise(dominant_genera = paste0(genus, collapse = "; "), .groups = "drop")

dmg_only_core_summary <- dmg_only_sample_summary %>%
  group_by(core) %>%
  summarise(
    n_samples = n_distinct(label),
    min_age_y_bp = min(y_bp, na.rm = TRUE),
    max_age_y_bp = max(y_bp, na.rm = TRUE),
    total_damaged_abundance = sum(total_damaged_abundance, na.rm = TRUE),
    median_n_genus_damaged = median(n_genus_damaged, na.rm = TRUE),
    median_n_species_damaged = median(n_species_damaged, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(dmg_only_dom_phylum, by = "core") %>%
  left_join(dmg_only_dom_genus, by = "core") %>%
  arrange(core)

write_tsv(dmg_only_core_summary, file.path(dmg_only_res_dir, "damaged_only_core_summary.tsv"))

# Figure 1: number of taxa through age (Damaged only)
dmg_only_n_taxa_long <- dmg_only_data %>%
  group_by(label) %>%
  summarise(
    species = n_distinct(species[!grepl("^Unclassified", species)]),
    genus = n_distinct(genus[!grepl("^Unclassified", genus)]),
    family = n_distinct(family[!grepl("^Unclassified", family)]),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(species, genus, family), names_to = "taxonomic_level", values_to = "n_taxa") %>%
  left_join(dmg_only_meta, by = "label") %>%
  mutate(
    core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")),
    taxonomic_level = factor(taxonomic_level, levels = c("family", "genus", "species"))
  )

dmg_only_n_taxa_plot <- ggplot(dmg_only_n_taxa_long, aes(x = y_bp, y = n_taxa, color = taxonomic_level)) +
  geom_point(size = 0.9, alpha = 0.9) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  facet_wrap(~core, scales = "free_x") +
  coord_flip() +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = "Age (y bp)", y = "Number of taxa (Damaged only)", color = "Taxonomic level") +
  theme(legend.position = "bottom")

png(file.path(dmg_only_plot_dir, "damaged_only_n_taxa.png"), width = 14, height = 8, units = "in", res = 600)
dmg_only_n_taxa_plot
dev.off()

# Figure 2: phylum relative abundance through age (Damaged only)
dmg_only_phylum_rel <- dmg_only_data %>%
  group_by(label, core, y_bp, phylum) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(label) %>%
  mutate(rel_abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(phylum) %>%
  summarise(mean_rel = mean(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(phylum_group = ifelse(mean_rel >= 0.02, phylum, "Other")) %>%
  select(phylum, phylum_group) %>%
  right_join(
    dmg_only_data %>%
      group_by(label, core, y_bp, phylum) %>%
      summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      group_by(label) %>%
      mutate(rel_abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
      ungroup(),
    by = "phylum"
  ) %>%
  mutate(phylum_group = ifelse(is.na(phylum_group), "Other", phylum_group)) %>%
  group_by(label, core, y_bp, phylum_group) %>%
  summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")))

dmg_only_phylum_rel_plot <- ggplot(dmg_only_phylum_rel, aes(x = y_bp, y = rel_abundance, fill = phylum_group)) +
  geom_area(position = "stack", color = "black", linewidth = 0.15) +
  facet_wrap(~core, scales = "free_x") +
  coord_flip() +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Age (y bp)", y = "Relative abundance within Damaged fraction", fill = "Phylum") +
  theme(legend.position = "bottom")

png(file.path(dmg_only_plot_dir, "damaged_only_phylum_rel_abundance.png"), width = 14, height = 8, units = "in", res = 600)
dmg_only_phylum_rel_plot
dev.off()

# Figure 3: heatmap of top genera across samples (Damaged only)
dmg_only_top_genera <- dmg_only_data %>%
  group_by(genus) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    n_samples_present = n_distinct(label[abundance > 0]),
    .groups = "drop"
  ) %>%
  mutate(prevalence_fraction = n_samples_present / ifelse(dmg_only_n_samples > 0, dmg_only_n_samples, NA_real_)) %>%
  arrange(desc(total_abundance), desc(prevalence_fraction)) %>%
  slice_head(n = 25)

dmg_only_heat_df <- dmg_only_data %>%
  filter(genus %in% dmg_only_top_genera$genus) %>%
  group_by(label, core, y_bp, genus) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(label) %>%
  mutate(rel_abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sample_order = paste(core, formatC(y_bp, width = 12, digits = 0, format = "f", flag = "0"), label, sep = "|"))

dmg_only_sample_levels <- dmg_only_meta %>%
  mutate(sample_order = paste(core, formatC(y_bp, width = 12, digits = 0, format = "f", flag = "0"), label, sep = "|")) %>%
  arrange(core, y_bp) %>%
  pull(sample_order)

dmg_only_heat_df <- dmg_only_heat_df %>%
  mutate(
    sample_order = factor(sample_order, levels = unique(dmg_only_sample_levels)),
    genus = fct_reorder(genus, abundance, .fun = sum, .desc = TRUE)
  )

dmg_only_top_genera_heatmap_plot <- ggplot(dmg_only_heat_df, aes(x = sample_order, y = genus, fill = log10(rel_abundance + 1e-6))) +
  geom_tile() +
  scale_fill_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#08306b", midpoint = -3, name = expression(log[10](Rel.~abundance + 10^{-6}))) +
  labs(x = "Samples ordered by core and age", y = "Genus") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

png(file.path(dmg_only_plot_dir, "damaged_only_top_genera_heatmap.png"), width = 15, height = 9, units = "in", res = 600)
dmg_only_top_genera_heatmap_plot
dev.off()

# Figure 4: top taxa barplot (Damaged only)
dmg_only_top_taxa_bar <- dmg_only_data %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(prop_total = total_abundance / sum(total_abundance, na.rm = TRUE)) %>%
  arrange(desc(prop_total)) %>%
  slice_head(n = 20) %>%
  mutate(genus = fct_reorder(genus, prop_total))

dmg_only_top_taxa_bar_plot <- ggplot(dmg_only_top_taxa_bar, aes(x = genus, y = prop_total)) +
  geom_col(fill = "#4C78A8", color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Genus", y = "Share of Damaged abundance", title = "Top genera in Damaged fraction")

png(file.path(dmg_only_plot_dir, "damaged_only_top_taxa_barplot.png"), width = 10, height = 7, units = "in", res = 600)
dmg_only_top_taxa_bar_plot
dev.off()

# Optional figure: domain share by core
dmg_only_domain_core <- dmg_only_data %>%
  group_by(core, domain) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(core) %>%
  mutate(rel_abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")))

dmg_only_domain_core_plot <- ggplot(dmg_only_domain_core, aes(x = core, y = rel_abundance, fill = domain)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Core", y = "Relative abundance (Damaged fraction)", fill = "Domain") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

png(file.path(dmg_only_plot_dir, "damaged_only_domain_share_by_core.png"), width = 8, height = 5.5, units = "in", res = 600)
dmg_only_domain_core_plot
dev.off()

# Conservative preliminary notes for manuscript drafting support
dmg_only_top_phyla_text <- dmg_only_data %>%
  group_by(phylum) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 5) %>%
  pull(phylum) %>%
  paste(collapse = ", ")

dmg_only_top_genera_text <- dmg_only_data %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 5) %>%
  pull(genus) %>%
  paste(collapse = ", ")

dmg_only_median_species <- median(dmg_only_sample_summary$n_species_damaged, na.rm = TRUE)
dmg_only_median_genus <- median(dmg_only_sample_summary$n_genus_damaged, na.rm = TRUE)

dmg_only_notes <- c(
  "- This section is a preliminary descriptive summary of the authenticated/damaged fraction only (is_dmg == 'Damaged').",
  "- Relative abundance values describe compositional share within the Damaged subset and do not imply absolute abundance in the original environment.",
  "- The Damaged fraction may reflect both original ecology and preservation/filtering processes; interpretations should remain conservative.",
  "- This section is intended to describe data structure and broad patterns before any network-based analyses.",
  paste0("- Across Damaged samples, median detected richness was ", round(dmg_only_median_genus, 1), " genera and ", round(dmg_only_median_species, 1), " species per sample."),
  paste0("- Most abundant Damaged phyla (overall): ", dmg_only_top_phyla_text, "."),
  paste0("- Most abundant Damaged genera (overall): ", dmg_only_top_genera_text, ".")
)

writeLines(dmg_only_notes, con = file.path(dmg_only_res_dir, "damaged_only_results_notes.txt"))
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
  filter(is_dmg == "Damaged") %>%
  mutate(phylum = paste(phylum, class, order, sep = "; ")) %>%
  group_by(is_dmg, core, domain, phylum) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  group_by(is_dmg, core) %>%
  mutate(abundance = abundance/sum(abundance)) %>%
  ungroup() %>%
  group_by(phylum) %>%
  mutate(total_abundance = mean(abundance)) %>%
  ungroup() %>%
  mutate(phylum = ifelse(total_abundance >= 0.005, phylum, "Other")) %>%
  group_by(is_dmg, core, domain, phylum) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  separate(phylum, into = c("phylum", "class", "order"), sep = "\\s*;\\s*", fill = "right", extra = "merge") %>%
  mutate(
    class = ifelse(is.na(class), "Other", class),
    order = ifelse(is.na(order), "", order)
  ) %>%
  mutate(class = paste(class, order, sep = "; ")) %>%
  mutate(class = ifelse(grepl("^Other", class), "Other", class)) %>%
  group_by(class) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  group_by(phylum) %>%
  mutate(total_abundance_p = sum(abundance)) %>%
  ungroup() %>%
  mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  mutate(class = { 
    tmp <- fct_reorder(class, total_abundance)
    if ("Other" %in% levels(tmp)) fct_relevel(tmp, "Other", after = 0) else tmp
  }) %>%
  mutate(phylum = { 
    tmp <- fct_reorder(phylum, total_abundance_p)
    if ("Other" %in% levels(tmp)) fct_relevel(tmp, "Other", after = 0) else tmp
  }) %>%
  mutate(phylum = fct_rev(phylum)) %>%
  ggplot(., aes(x = class, y = abundance, fill = phylum))+
    geom_bar(position = "dodge", stat = "identity", color = "black", linewidth = 0.2)+
    facet_nested(domain~core, scales = "free_y", space = "free_y")+
    coord_flip()+
    # scale_y_sqrt()+
    guides(fill = guide_legend(ncol = 1))+
    scale_fill_manual(values= c(paired_genus[1:21], "grey"))+
    labs(x = "Tax (class and order)", y = "Relative abundance", fill = "Phylum", size = "Number of reads")

png(file = "plots/microbial/preliminary/order_damage_toptax.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()




# Correlaciones entre "orders" (misma lógica que para species, aplicado a niveles 'order')
corr_wide_orders <- tax_data %>%
  filter(is_dmg == "Damaged") %>%
  group_by(is_dmg, core, label, y_bp, domain, phylum, class, order) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  # filter(core == "GeoB25202_R1") %>%
  # tratar NA en order y calcular abundancia total por order
  mutate(order = ifelse(is.na(order) | order == "", "Other", order)) %>%
  group_by(order) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  # quedarnos con orders con suficiente señal (ajusta el umbral si hace falta)
  filter(total_abundance > 10000) %>%

  # calcular abundancias relativas por muestra (label) entre los orders seleccionados
  group_by(label) %>%
  mutate(rel_ab = abundance / sum(abundance)) %>%
  ungroup() %>%
  select(label, order, rel_ab) %>%
  pivot_wider(names_from = order, values_from = rel_ab, values_fill = 0)

# convertir a matriz (filas = muestras, columnas = orders)
corr_mat_input_o <- corr_wide_orders %>% as.data.frame()
rownames(corr_mat_input_o) <- corr_mat_input_o$label
mat_o <- as.matrix(corr_mat_input_o[ , setdiff(colnames(corr_mat_input_o), "label")])

# matriz de correlaciones (Spearman)
corr_mat_o <- cor(mat_o, use = "pairwise.complete.obs", method = "spearman")

# calcular matriz de p-valores (pairwise Spearman cor.test)
ords <- colnames(mat_o)
n_ord <- length(ords)
p_mat_o <- matrix(NA_real_, nrow = n_ord, ncol = n_ord, dimnames = list(ords, ords))
for (i in seq_len(n_ord)) {
  for (j in i:n_ord) {
    xi <- mat_o[, i]
    yj <- mat_o[, j]
    idx <- complete.cases(xi, yj)
    if (sum(idx) < 3) {
      pval <- NA_real_
    } else {
      pval <- tryCatch(cor.test(xi[idx], yj[idx], method = "spearman")$p.value,
                       error = function(e) NA_real_)
    }
    p_mat_o[i, j] <- pval
    p_mat_o[j, i] <- pval
  }
}
diag(p_mat_o) <- 0

# ordenar por clustering jerárquico para la visualización (manejar posible error si hay NA)
ord_hc <- tryCatch({
  hc_o <- hclust(as.dist(1 - corr_mat_o))
  hc_o$order
}, error = function(e) {
  seq_len(ncol(corr_mat_o))
})
orders_order <- colnames(corr_mat_o)[ord_hc]

# pasar a formato largo y aplicar orden, añadiendo p-valores
corr_long_o <- as.data.frame(corr_mat_o) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "corr")

p_long_o <- as.data.frame(p_mat_o) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "pval")

corr_long_o <- corr_long_o %>%
  left_join(p_long_o, by = c("Var1", "Var2")) %>%
  mutate(
    Var1 = factor(Var1, levels = orders_order),
    Var2 = factor(Var2, levels = orders_order),
    # conservar solo correlaciones significativas (p < 0.05)
    corr_sig = ifelse(!is.na(pval) & pval < 0.05, corr, NA_real_)
  )

# dibujar y guardar heatmap: solo tiles con correlación significativa estarán coloreadas
# png("plots/microbial/preliminary/nitrososphaerales_orders_corr_heatmap.png", width = 8, height = 8, units = "in", res = 300)
ggplot(corr_long_o, aes(x = Var2, y = Var1, fill = corr_sig)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1), na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "", y = "", fill = "Spearman\ncorr", title = "Correlación entre orders (GeoB25202_R1, Damaged)") 
  # geom_text(data = subset(corr_long_o, !is.na(corr_sig)), aes(label = round(corr, 2)), size = 3, color = "black")
# dev.off()
     

     
     



