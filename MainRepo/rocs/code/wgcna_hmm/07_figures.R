#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggh4x)



  library(forcats)


})

# -----------------------------------------------------------------------------
# Exploratory WGCNA figure objects
# -----------------------------------------------------------------------------

source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")
setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")

cdata <- read_tsv("data/metadata_v5.tsv")






main_states <- read_tsv("results/wgcna_hmm_leiden/main/hmm_states_main.tsv")

states_def <- main_states %>%
  filter(core == "GeoB25202_R1") %>%
  group_by(state) %>%
  summarise(mean_sst = mean(sst, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_sst) %>%
  mutate(state_definition = c("state_01", "state_02", "state_03", "state_04", "state_05")) %>%
  select(state, state_definition)

main_states <- main_states %>%
  select(-d18O, -sst, -label) %>%
  left_join(states_def, by = "state")

color_states <- c("#3f74b3", "#59b3a9", "#fdcd7b", "#f16943", "#9d0241")
names(color_states) <- c("state_01", "state_02", "state_03", "state_04", "state_05")


# 1) Construir intervalos por muestra (midpoints) dentro de cada core
states_rect <- main_states %>%
  mutate(y_bp = age_kyr * 1000) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  arrange(core, y_bp) %>%
  group_by(core) %>%
  mutate(
    x_left_raw  = (y_bp + lag(y_bp)) / 2,
    x_right     = (y_bp + lead(y_bp, default = last(y_bp))) / 2,
    # para el primer punto del core, arrancar en 0
    x_left      = if_else(row_number() == 1, 0, x_left_raw)
  ) %>%
  select(-x_left_raw) %>%
  ungroup()

# 2) Colapsar intervalos consecutivos con el mismo estado (para que no haya 1 rectángulo por muestra)
states_runs <- states_rect %>%
  arrange(core, y_bp) %>%
  group_by(core) %>%
  mutate(run_id = cumsum(state_definition != lag(state_definition, default = first(state_definition)))) %>%
  ungroup() %>%
  group_by(core, run_id, state_definition) %>%
  summarise(
    xmin = min(x_left),
    xmax = max(x_right),
    .groups = "drop"
  )

# 3) Plot: rectángulos cubriendo TODO el eje
plotA_full <- ggplot(states_runs, aes(fill = state_definition)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -1, ymax = 1),
            color = "black", linewidth = 0.15) +
  facet_nested(~core, scales = "free") +
  coord_flip() +
  scale_fill_manual(values = color_states) +
  scale_x_reverse(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = seq(0, 6e06, 1e4)
  ) +
  labs(x = "Age (y bp)", y = NULL, fill = "State") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"

  )

plotB <- cdata %>%
  filter(label %in% main_states$sample) %>%
  filter(grepl("GeoB25202_R1", core)) %>%
  ggplot(., aes(x = y_bp, y = sst)) +
    geom_point()+
    geom_line()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4)) +
    coord_flip()


png("states.png", width = 3000, height = 2000, res = 150)
egg::ggarrange(plotA_full, plotB, nrow = 1, widths = c(3, 1))
dev.off()






extended_states <- read_tsv("results/wgcna_hmm_leiden/extended/hmm_states_extended_selected.tsv") %>%
  select(sample, state) %>%
  left_join(states_def, by = "state") %>%
  left_join(cdata, by = c("sample" = "label"))

# extended_states %>%
#   mutate(y_bp = age_kyr * 1000) %>%
#   mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#   ggplot(., aes(x = y_bp, y = 1)) +
#     geom_point(aes(colour = as.factor(state))) +
#     geom_segment(aes(y = 1, yend = -1, colour = as.factor(state)), size = 2.5, lineend = "butt") +
#     facet_nested(~core, scales = "free")+
#     coord_flip() +
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))


# 1) Construir intervalos por muestra (midpoints) dentro de cada core
states_rect <- extended_states %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  arrange(core, y_bp) %>%
  group_by(core) %>%
  mutate(
    x_left_raw  = (y_bp + lag(y_bp)) / 2,
    x_right     = (y_bp + lead(y_bp, default = last(y_bp))) / 2,
    # para el primer punto del core, arrancar en 0
    x_left      = if_else(row_number() == 1, 0, x_left_raw)
  ) %>%
  select(-x_left_raw) %>%
  ungroup()

# 2) Colapsar intervalos consecutivos con el mismo estado (para que no haya 1 rectángulo por muestra)
states_runs <- states_rect %>%
  arrange(core, y_bp) %>%
  group_by(core) %>%
  mutate(run_id = cumsum(state_definition != lag(state_definition, default = first(state_definition)))) %>%
  ungroup() %>%
  group_by(core, run_id, state_definition) %>%
  summarise(
    xmin = min(x_left),
    xmax = max(x_right),
    .groups = "drop"
  ) %>%
  mutate(transparency = ifelse(core %in% c("ST8", "ST13", "GeoB25202_R1") & xmax <= 150000, 1, 0))

# 3) Plot: rectángulos cubriendo TODO el eje
plotA_full <- ggplot(states_runs, aes(fill = state_definition)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -1, ymax = 1),
            color = "black", linewidth = 0.15) +
  # geom_rect(aes(xmin = 150000, xmax = 600000, ymin = -1, ymax = 1),
  #           fill = "white", linewidth = 0.15, alpha = 0.02) +
  facet_nested(~core, scales = "free") +
  coord_flip() +
  scale_fill_manual(values = color_states) +
  scale_x_reverse(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = seq(0, 6e06, 1e4)
  ) +
  labs(x = "Age (y bp)", y = NULL, fill = "State") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

plotB <- cdata %>%
  filter(label %in% extended_states$sample) %>%
  filter(grepl("GeoB25202_R1", core)) %>%
  ggplot(., aes(x = y_bp, y = sst)) +
    geom_point()+
    geom_line()+
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4)) +
    coord_flip()+
    labs(x = NULL, y = "SST (°C)") +
    theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

plotC <- cdata %>%
  filter(label %in% extended_states$sample) %>%
  filter(grepl("GeoB25202_R1", core)) %>%
  ggplot(., aes(x = y_bp, y = -mis)) +
    geom_point()+
    geom_line()+
    labs(x = NULL, y = "d18O (‰)") +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4)) +
    coord_flip()+
    theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

png("extended_states.png", width = 3000, height = 3000, res = 150)
egg::ggarrange(plotA_full, plotB, nrow = 1, widths = c(5, 1))
dev.off()






## Eigengenes
training_eigengenes <- read_tsv("results/wgcna_hmm_leiden/main/module_eigengenes_training.tsv")

training_eigengenes %>% 
  pivot_longer(cols = starts_with("ME"), names_to = "module", values_to = "eigengene") %>%
  left_join(cdata, by = c("sample" = "label")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = eigengene)) +
    geom_point(aes(colour = core)) +
    geom_line(aes(colour = core)) +
    geom_smooth(method = "loess", se = TRUE, span = 0.1, color = "black") +
    facet_nested(~module, scales = "free")+
    coord_flip() +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))
  
training_eigengenes %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module", values_to = "eigengene") %>%
  left_join(cdata, by = c("sample" = "label")) %>%
  # filter(grepl("GeoB25202", core)) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = sst, y = eigengene))+
    geom_point()+
    facet_wrap(~module, scales = "free")+
    geom_smooth(method = "lm", se = TRUE, color = "black")

eigen_sst <- training_eigengenes %>%
  left_join(cdata, by = c("sample" = "label")) %>%
  filter(grepl("GeoB25202", core))

cor.test(eigen_sst$sst, eigen_sst$MEbrown, method = "spearman", use = "complete.obs")
cor.test(eigen_sst$sst, eigen_sst$MEblue, method = "spearman", use = "complete.obs")
cor.test(eigen_sst$sst, eigen_sst$MEgreen, method = "spearman", use = "complete.obs")
cor.test(eigen_sst$sst, eigen_sst$MEgrey, method = "spearman", use = "complete.obs")
cor.test(eigen_sst$sst, eigen_sst$MEturquoise, method = "spearman", use = "complete.obs")
cor.test(eigen_sst$sst, eigen_sst$MEyellow, method = "spearman", use = "complete.obs")


projection_eigengenes <- read_tsv("results/wgcna_hmm_leiden/extended/module_eigengenes_extended_projection.tsv")

plotA <- projection_eigengenes %>% 
  # mutate(y_bp = age_kyr * 1000) %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module", values_to = "eigengene") %>%
  # filter(!grepl("grey", module)) %>%
  left_join(cdata, by = c("sample" = "label")) %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = eigengene)) +
    geom_point(aes(colour = core)) +
    geom_line(aes(colour = core)) +
    geom_smooth(method = "loess", se = TRUE, span = 0.05, color = "black") +
    facet_nested(~module*core, scales = "free")+
    coord_flip() +
    theme(legend.position = "left") +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))

png("extended_eigengenes.png", width = 3000, height = 3000, res = 150)
egg::ggarrange(plotA, plotB, nrow = 1, widths = c(5, 1))
dev.off()



modules_xrf <- read_tsv("results/wgcna_hmm_leiden/main/loco_stable_associations.tsv")

df_hm <- modules_xrf %>%
  mutate(
    sig = case_when(
      p_full < 0.001 ~ "***",
      p_full < 0.01  ~ "**",
      p_full < 0.05  ~ "*",
      TRUE ~ ""
    ),
    module = factor(module),
    proxy  = factor(proxy)
  )

p <- ggplot(df_hm, aes(x = proxy, y = module, fill = beta_full)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Proxy", y = "Module", fill = "Correlation (beta)")
p






data <- read_tsv("results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz")
module_assignment <- read_tsv("results/wgcna_hmm_leiden/main/module_assignments.tsv")


plotA <- data %>%
  filter(is_dmg == "Damaged" & core != "ST5") %>%
  # filter(!phylum %in% c("p__Asgardarchaeota", "p__Aerophobota")) %>%
  select(label, abundance, domain, phylum, class, order, family, genus, species, subspecies) %>%
  inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
  mutate(tax = paste(domain, phylum, sep = ";")) %>%
  group_by(label, tax, module) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(label, tax, module, fill = list(abundance = 0)) %>%
  left_join(cdata, by = c("label" = "label")) %>%
  # filter(module %in% c("blue", "green")) %>%
  group_by(tax) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(tax = fct_reorder(tax, desc(total_abundance))) %>%
  # filter(grepl("GeoB25202_R1", core)) %>%
  # group_by(module, label, core, y_bp, domain, phylum, class, order, is_dmg) %>%
  # summarise(abundance = sum(abundance), .groups = "drop") %>%
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  ggplot(., aes(x = y_bp, y = abundance, fill = tax)) +
    geom_area(stat = "identity", position = "stack") +
    facet_nested(~module*core, scales = "free") +
    # scale_y_sqrt()+
    scale_fill_manual(values = paired_genus)+
    guides(fill = guide_legend(title = "Taxonomic group", ncol = 1)) +
    coord_flip() +
    theme(legend.position = "left") +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))


png("extended_abundances.png", width = 4000, height = 2000, res = 150)
egg::ggarrange(plotA, plotB, nrow = 1, widths = c(5, 1))
dev.off()


plotA <- {
  # compute abundance per tax within the blue module only
  tax_blue <- data %>%
    filter(is_dmg == "Damaged" & core != "ST5") %>%
    inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
    filter(module == "leiden1") %>%
    mutate(tax = paste(domain, phylum, class, order, family, genus, sep = ";")) %>%
    group_by(tax) %>%
    summarise(blue_abundance = sum(abundance), .groups = "drop")

  # main plot pipeline, joining the blue-specific abundance and ordering by it
  data %>%
    filter(is_dmg == "Damaged" & core != "ST5") %>%
    select(abundance, domain, phylum, class, order, family, genus, species, subspecies) %>%
    inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
    filter(!module %in% c("grey")) %>%
    mutate(tax = paste(domain, phylum, class, order, family, genus, sep = ";")) %>%
    group_by(tax, module) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    left_join(tax_blue, by = "tax") %>%
    mutate(blue_abundance = replace_na(blue_abundance, 0)) %>%
    mutate(tax = fct_reorder(tax, blue_abundance)) %>%
    ggplot(., aes(x = tax, y = abundance)) +
      geom_point() +
      facet_nested(~module, scales = "free") +
      guides(fill = guide_legend(title = "Taxonomic group", ncol = 1)) +
      coord_flip()
}

png("plot_test.png", width = 3000, height = 2000, res = 150)
plotA
dev.off()


plotA <- data %>%
    filter(is_dmg == "Damaged" & core != "ST5") %>%
    inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
    # filter(module == "leiden6") %>%
    mutate(tax = paste(domain, phylum, class, order, family, genus, sep = ";")) %>%
    group_by(tax, module) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    group_by(module) %>%
    mutate(abundance = abundance / sum(abundance)) %>%
    ungroup() %>%
    group_by(tax) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(tax = fct_reorder(tax, desc(total_abundance))) %>%
    filter(grepl("Colwellia", tax)) %>% 
    ggplot(., aes(x = module, y = abundance, fill = tax))+
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = paired_genus)+
      guides(fill = guide_legend(title = "Taxonomic group", ncol = 2))

png("plot_test.png", width = 3000, height = 2000, res = 150)
plotA
dev.off()

data %>%
  filter(grepl("Polar", genus)) %>%
  select(domain, phylum, class, order, family, genus, species, subspecies) %>%
  head()


data %>%
  filter(is_dmg == "Damaged" & core != "ST5") %>%
  select(abundance, domain, phylum, class, order, family, genus, species, subspecies) %>%
  inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
  filter(module != "grey") %>%
  filter(module %in% c("blue", "green")) %>%
  mutate(tax = paste(phylum, class, sep = ";")) %>%
  group_by(domain, tax, module) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  group_by(domain, module) %>%
  mutate(prop_abundance = abundance / sum(abundance)) %>%
  ungroup() %>%
  group_by(tax) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(tax = fct_reorder(tax, desc(total_abundance))) %>%
  ggplot(., aes(x = module, y = prop_abundance, fill = tax)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = paired_genus)+
    facet_nested(~domain, scales = "free") +
    guides(fill = guide_legend(title = "Taxonomic group", ncol = 2)) 
    # theme(legend.position = "none")
  
  







plotB <- data %>% 
  filter(is_dmg == "Damaged") %>%
  # select(label, abundance, domain, phylum, class, order, family, genus, species, subspecies, is_dmg, core, y_bp) %>%
  inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
  filter(module == "brown") %>% 
  filter(grepl("GeoB25202_R1", core)) %>%
  # group_by(module, label, core, y_bp, domain, phylum, class, order, is_dmg) %>%
  # summarise(abundance = sum(abundance), .groups = "drop") %>%
  ggplot(., aes(x = y_bp, y = temp)) +
    geom_line() +
    facet_nested(~core, scales = "free") +
    scale_y_sqrt()+
    coord_flip() +
    f
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))

plotC <- data %>% 
  filter(is_dmg == "Damaged") %>%
  # select(label, abundance, domain, phylum, class, order, family, genus, species, subspecies, is_dmg, core, y_bp) %>%
  inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
  filter(module == "brown") %>% 
  filter(grepl("GeoB25202_R1", core)) %>%
  # group_by(module, label, core, y_bp, domain, phylum, class, order, is_dmg) %>%
  # summarise(abundance = sum(abundance), .groups = "drop") %>%
  ggplot(., aes(x = y_bp, y = initial)) +
    geom_line() +
    facet_nested(~core, scales = "free") +
    scale_y_sqrt()+
    coord_flip() +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))

egg::ggarrange(plotC, plotA, plotB, nrow = 1, widths = c(1, 2, 1))


data %>%
  filter(is_dmg == "Damaged") %>%
  select(label, abundance, domain, phylum, class, order, family, genus, species, subspecies, is_dmg, core, y_bp) %>%
  inner_join(module_assignment, by = c("subspecies" = "taxon")) %>%
  group_by(module, domain, phylum, class, order, family, genus, species, subspecies, is_dmg) %>%
  summarise(abundance = sum(abundance), n_labels = n_distinct(label), .groups = "drop") %>%
  filter(module == "brown") %>%
  arrange(desc(abundance)) %>%
  View()



# fingerprint <- read_tsv("results/wgcna_hmm/main/state_fingerprints_main.tsv")



















cdataxrf <- read_tsv("results/wgcna_hmm/main/module_xrf_dataset_training.tsv")
cdataxrf %>% 
  mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
  left_join(extended_states %>% select(label = sample, state), by = "label") %>%
  ggplot(., aes(x = y_bp))+
    geom_segment(aes(y = 0, yend = Inf, xend = y_bp, color = as.factor(state))) +
    geom_line(aes(y = ca/ti), color = "black") +
    geom_line(aes(y = ba/ti), color = "red") +
    facet_nested(~core, scales = "free")+
    coord_flip() +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 1e4))
  














paths <- list(
  sample_metadata = file.path(BASE, "results/microbial/wgcna/sample_metadata_stage1.tsv"),
  harmonized_eigengenes = file.path(BASE, "results/microbial/wgcna/harmonized/wgcna/module_eigengenes.tsv"),
  harmonized_hmm_states = file.path(BASE, "results/microbial/wgcna/harmonized/hmm/hmm_states.tsv"),
  harmonized_state_fingerprints = file.path(BASE, "results/microbial/wgcna/harmonized/hmm/state_fingerprints.tsv"),
  harmonized_module_xrf_dataset = file.path(BASE, "results/microbial/wgcna/harmonized/downstream_training_only/module_xrf_dataset.tsv"),
  harmonized_module_assignments = file.path(BASE, "results/microbial/wgcna/harmonized/wgcna/module_assignments.tsv"),
  prok_taxa_metadata = file.path(BASE, "results/microbial/wgcna/prokaryotes_taxa_metadata.tsv")
)

sample_meta <- fread(paths$sample_metadata)
harmonized_me <- fread(paths$harmonized_eigengenes)
harmonized_hmm <- fread(paths$harmonized_hmm_states)
state_fingerprints <- fread(paths$harmonized_state_fingerprints)
module_xrf_dataset <- fread(paths$harmonized_module_xrf_dataset)
module_assignments <- fread(paths$harmonized_module_assignments)
prok_taxa_metadata <- fread(paths$prok_taxa_metadata)

core_levels <- c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")




sample_meta %>%
  left_join(harmonized_me, by = c("label" = "sample")) %>%
  group_by(core) %>%
  summarise(
    across(
      starts_with("ME"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        n = ~ sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) 


sample_meta %>%
  left_join(harmonized_me, by = c("label" = "sample")) %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module", values_to = "eigengene") %>%
  ggplot(., aes(x = y_bp, y = eigengene, colour = core)) +
  geom_point() +
  facet_wrap(~module, scales = "free_y") +
  labs(
    x = "Age (kyr BP)",
    y = "Eigengene",
    colour = "Core",
    title = "Harmonized module eigengenes through time"
  ) +
  coord_flip() +
  scale_x_reverse() +
  theme_bw(base_size = 10)

harmonized_me %>% head()




dt <- sample_meta %>%
  left_join(harmonized_me, by = c("label" = "sample"))



# log_msg("Residualising eigengenes within core...")
me_cols <- grep("^ME(?!grey)", names(new_MEs), value = TRUE, perl = TRUE)
# log_msg(sprintf("  Using %d module eigengenes: %s", length(me_cols),
                # paste(me_cols, collapse = ", ")))


for (col in me_cols) {
  core_means <- dt[, .(core_mean = mean(get(col), na.rm = TRUE)), by = core]
  dt[core_means, paste0(col, "_resid") := get(col) - i.core_mean, on = "core"]
}
resid_cols <- paste0(me_cols, "_resid")


dt %>%
  select(label, core, y_bp, all_of(resid_cols)) %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module", values_to = "eigengene") %>%
  ggplot(., aes(x = y_bp, y = eigengene, colour = core)) +
  geom_point() +
  facet_wrap(~module, scales = "free_y") +
  labs(
    x = "Age (kyr BP)",
    y = "Eigengene",
    colour = "Core",
    title = "Harmonized module eigengenes through time"
  ) +
  coord_flip() +
  scale_x_reverse() +
  theme_bw(base_size = 10)





# 1) Eigengene residuals through time ------------------------------------------

me_cols <- grep("^ME", names(harmonized_me), value = TRUE)

harmonized_me_long <- harmonized_me %>%
  inner_join(
    sample_meta %>%
      transmute(sample = label, core, y_bp, age_kyr = y_bp / 1000),
    by = "sample"
  ) %>%
  pivot_longer(cols = all_of(me_cols), names_to = "module", values_to = "eigengene") %>%
  mutate(
    module = str_remove(module, "^ME"),
    core = factor(core, levels = core_levels)
  ) %>%
  group_by(core, module) %>%
  mutate(eigengene_resid = eigengene - mean(eigengene, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(eigengene_resid), !is.na(age_kyr), !is.na(core))

plot_eigengene_residuals_time <- harmonized_me_long %>%
  filter(grepl("green", module))  %>%
  ggplot(aes(x = age_kyr, y = eigengene_resid, colour = core)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey35") +
    geom_point(alpha = 0.45, size = 0.7) +
    geom_smooth(
      aes(group = 1),
      method = "loess",
      se = FALSE,
      colour = "black",
      linewidth = 0.45
    ) +
    geom_smooth(
      method = "loess",
      se = FALSE,
      linewidth = 0.35,
      alpha = 0.7
    ) +
    coord_flip()+
    facet_nested(~module*core, scales = "free") +
    scale_x_reverse() +
    labs(
      x = "Age (kyr BP)",
      y = "Eigengene residual (within-core centered)",
      colour = "Core",
      title = "Harmonized module eigengene residuals through time"
    ) +
    theme_bw(base_size = 10)

# 2) Correlations between module eigengenes and XRF ---------------------------

xrf_candidate_cols <- names(module_xrf_dataset)[
  names(module_xrf_dataset) %in% c(
    "ba", "ti", "ca", "zr", "rb", "sr", "al", "fe", "mn", "br", "p", "si",
    "ba_ti_ratio", "ca_ti_ratio", "ti_ca_ratio", "zr_rb_ratio", "zr_sr_ratio",
    "ti_al_ratio", "fe_al_ratio", "mn_fe_ratio", "br_ti_ratio", "p_ti_ratio", "si_ti_ratio"
  )
]

me_xrf_cols <- grep("^ME", names(module_xrf_dataset), value = TRUE)

module_xrf_cor <- expand_grid(
  module = me_xrf_cols,
  xrf = xrf_candidate_cols
) %>%
  rowwise() %>%
  mutate(
    cor_value = cor(module_xrf_dataset[[module]], module_xrf_dataset[[xrf]], use = "pairwise.complete.obs"),
    n_complete = sum(is.finite(module_xrf_dataset[[module]]) & is.finite(module_xrf_dataset[[xrf]]))
  ) %>%
  ungroup() %>%
  mutate(
    module = str_remove(module, "^ME"),
    xrf = factor(xrf, levels = unique(xrf))
  )

plot_module_xrf_correlations <- module_xrf_cor %>%
  ggplot(aes(x = xrf, y = fct_rev(module), fill = cor_value)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    x = "XRF variable",
    y = "Module",
    fill = "Pearson r",
    title = "Module eigengene–XRF correlations (harmonized training-only dataset)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# 3) Classification of samples into states across all cores --------------------

# Use harmonized HMM state calls (contains all cores + labels).
hmm_timeline <- harmonized_hmm %>%
  transmute(
    sample,
    core = factor(core, levels = core_levels),
    age_kyr,
    state = factor(label)
  ) %>%
  filter(!is.na(core), !is.na(age_kyr), !is.na(state)) %>%
  arrange(core, age_kyr) %>%
  group_by(core) %>%
  mutate(age_kyr_next = lead(age_kyr)) %>%
  ungroup()

plot_states_all_cores_timeline <- hmm_timeline %>%
  ggplot(aes(y = core, colour = state)) +
  geom_segment(
    data = ~ filter(.x, !is.na(age_kyr_next)),
    aes(x = age_kyr, xend = age_kyr_next, yend = core),
    linewidth = 2.7,
    lineend = "butt",
    alpha = 0.9
  ) +
  geom_point(aes(x = age_kyr), size = 1.1, alpha = 0.95) +
  scale_x_reverse() +
  coord_flip() +
  labs(
    x = "Age (kyr BP)",
    y = "Core",
    colour = "State",
    title = "State classification timeline across all cores"
  ) +
  theme_bw(base_size = 10)

# 4) Taxonomic composition of each module --------------------------------------

module_taxonomy <- module_assignments %>%
  filter(!is.na(module), module != "grey") %>%
  left_join(
    prok_taxa_metadata %>% select(taxon, domain, phylum),
    by = "taxon"
  ) %>%
  mutate(
    domain = replace_na(domain, "Unknown"),
    phylum = replace_na(phylum, "Unknown")
  ) %>%
  count(module, phylum, name = "n_taxa") %>%
  group_by(module) %>%
  mutate(prop_taxa = n_taxa / sum(n_taxa)) %>%
  ungroup()

top_phyla <- module_taxonomy %>%
  group_by(phylum) %>%
  summarise(total = sum(n_taxa), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 12) %>%
  pull(phylum)

module_taxonomy_plot_df <- module_taxonomy %>%
  mutate(phylum_plot = if_else(phylum %in% top_phyla, phylum, "Other")) %>%
  group_by(module, phylum_plot) %>%
  summarise(prop_taxa = sum(prop_taxa), .groups = "drop")

plot_module_taxonomic_composition <- module_taxonomy_plot_df %>%
  ggplot(aes(x = module, y = prop_taxa, fill = phylum_plot)) +
  geom_col(width = 0.85) +
  labs(
    x = "Module",
    y = "Proportion of assigned taxa",
    fill = "Phylum",
    title = "Taxonomic composition of modules (non-grey taxa assignments)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# 5) Fingerprint of each state relative to modules -----------------------------

fingerprint_long <- state_fingerprints %>%
  pivot_longer(
    cols = starts_with("ME"),
    names_to = "module",
    values_to = "fingerprint"
  ) %>%
  mutate(module = str_remove(module, "^ME"))

plot_state_module_fingerprint <- fingerprint_long %>%
  ggplot(aes(x = module, y = fct_reorder(label, as.numeric(state), .fun = mean), fill = fingerprint)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    x = "Module",
    y = "State",
    fill = "Fingerprint",
    title = "State fingerprints across module eigengenes (harmonized HMM)"
  ) +
  theme_bw(base_size = 10)
