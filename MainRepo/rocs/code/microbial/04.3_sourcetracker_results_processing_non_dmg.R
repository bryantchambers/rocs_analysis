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



############################################
##### Process sourcetracker outputs  #######
############################################

#### Output sourcetracker results
cdata <- read.table(file = "data/metadata_v5.tsv", sep = "\t", header = T)

source_cdata <- read_tsv("data/sourcetracker/med-biomes-download.txt") %>%
    mutate(biome_subclass = gsub(":|-| +", "_", biome_subclass)) %>%
    mutate(biome_subclass = tolower(biome_subclass)) %>%
    mutate(biome_subclass = gsub(" +", "_", biome_subclass))

header_line <- readLines("results/microbial/sourcetracker/st_sp_median_wviruses-non-dmg/mixing_proportions.txt", n = 1) %>% sub("^#", "", .)
stout <- read.table(file = "results/microbial/sourcetracker/st_sp_median_wviruses-non-dmg/mixing_proportions.txt", sep = "\t", header = FALSE, col.names = unlist(strsplit(header_line, "\t"))) %>% 
    janitor::clean_names() %>%
    rename(label = sample_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -label, names_to = "source", values_to = "proportion") %>%
    left_join(source_cdata %>% select(biome_class, biome_subclass) %>% distinct(), by = c("source"="biome_subclass")) %>%
    mutate(biome_class = ifelse(source == "unknown", "Unknown", biome_class)) %>%
    mutate(
        source = gsub("^root_environmental_aquatic_", "", source),
        source = gsub("^root_environmental_terrestrial_soil_", "", source),
        source = gsub("^marine_", "", source),
        source = source %>%
            str_replace_all("_", " ") %>%
            str_to_lower() %>%
            str_to_title()) %>%
    mutate(biome_class = gsub("^root:", "", biome_class))


# Taxonomic profile at family level
sources <- stout %>%
    group_by(source) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    arrange(desc(proportion)) %>%
    pull(source) %>% unique()
paired_genus_add <- c(c("lightgrey", "#154360", "#FF5733", "#FFC300", "#1ABC9C"), paired_genus)
colors <- paired_genus_add[1:length(sources)]
names(colors) <- sources

plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source, biome_class, core)  %>%
    summarise(mean_proportion = mean(proportion), sd = sd(proportion), .groups = "drop") %>%
    group_by(source) %>%
    mutate(total_proportion = sum(mean_proportion)) %>%
    ungroup() %>%
    mutate(source = fct_reorder(source, total_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(sd = ifelse(sd > mean_proportion, mean_proportion-0.0000001, sd)) %>%
    mutate(biome_class = fct_relevel(biome_class, "Unknown", "Environmental:Aquatic:Marine", "Environmental:Aquatic:Freshwater", "Environmental:Terrestrial:Soil", "Host-associated:Plants")) %>%
    # filter(biome_class %in% c("Unknown", "Environmental:Aquatic:Marine")) %>%
    ggplot(., aes(x = source, y = mean_proportion, fill = source)) +
            geom_bar(stat = "identity") +
            geom_errorbar(aes(ymin = mean_proportion - sd, ymax = mean_proportion + sd), width = 0.2) +
            coord_flip() +
            scale_fill_manual(values = colors)+
            facet_nested(biome_class~core, scales = "free_y", space = "free_y") +
            # scale_y_sqrt(labels = scales::percent_format(), breaks = c(0, .05, .25, .50, .75))+
            scale_y_sqrt(labels = scales::percent_format())+
            labs(x = "", y = "Mean source proportion", fill = "Source")+
            # labs(x = "", y = "Mean source proportion", fill = "Source", title = "Mean source proportion by biome class and core")+
            theme(legend.position = "none")

# png(file = "plots/microbial/sourcetracker/01_proportion_by_core.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()


plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    mutate(source = fct_reorder(source, mean_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
        geom_area(stat = "identity", position = "stack") +
        labs(x = "Age (y bp)", y = "Source proportion", fill = "Source") +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4)) +
        facet_nested(~core, scales = "free", space = "fixed") +
        scale_fill_manual(values = colors) +
        coord_flip()+
        guides(fill = guide_legend(reverse = TRUE))

# png(file = "plots/microbial/sourcetracker/02_proportion_sample.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()


plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    filter(source %in% c("Polar", "Atlantic Ocean", "Estuary", "Pacific Ocean", "Mediterranean Sea")) %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "Age (y bp)", y = "Source proportion", fill = "Source")+
            scale_y_continuous(labels = scales::percent_format(), breaks = breaks_pretty(n = 3))+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
            facet_nested(.~source*core, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()+
            guides(fill = guide_legend(ncol = 1))
            
# png(file = "plots/microbial/sourcetracker/03_proportion_sample_facet.png", width = 30, height = 10, units = "in", res = 600)
plotA
# dev.off()


plotA <- stout %>%
    left_join(cdata) %>%
    filter(grepl("Polar|Atlantic Ocean", source)) %>%
    select(label, y_bp, source, proportion) %>%
    pivot_wider(names_from = source, values_from = proportion, values_fill = 0) %>%
    clean_names() %>%
    mutate(total = atlantic_ocean + polar) %>%
    mutate(ratio = log2((atlantic_ocean + 1e-10) / (polar + 1e-10))) %>%
    filter(!is.nan(ratio)) %>%
    left_join(cdata) %>%
    filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(
        ratio_pos = ifelse(ratio > 0, ratio, 0),
        ratio_neg = ifelse(ratio < 0, ratio, 0)
    ) %>%
    ggplot(aes(x = y_bp)) +
        geom_area(aes(y = ratio_neg, fill = "Polar-dominated")) +
        geom_area(aes(y = ratio_pos, fill = "Atlantic-dominated")) +
        geom_line(aes(y = ratio), color = "black", alpha = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        coord_flip() +
        facet_nested(~core, scales = "free", space = "fixed") +
        scale_x_reverse(
            labels = scales::label_number(scale_cut = scales::cut_short_scale()),
            breaks = seq(0, 5e06, 2e4)) +
        scale_fill_manual(
            name = NULL,
            breaks = c("Polar-dominated", "Atlantic-dominated"),
            values = c("Atlantic-dominated" = "#FF5733", "Polar-dominated" = "#154360")) +
        labs(x = "Age (y bp)", y = expression(log[2](Atlantic / Polar~source~proportion)))

# png(file = "plots/microbial/sourcetracker/04_atlantic_polar_ratio.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()





stout_wide <- stout %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    select(label, source, proportion) %>%
    pivot_wider(
    names_from = source,
    values_from = proportion,
    values_fill = 0
    )

# Calcular matriz de correlación y matriz de p-values usando cor.test para cada par de columnas
sources <- stout_wide %>% select(-label)
n <- ncol(sources)
cor_mat <- matrix(NA, n, n)
p_mat <- matrix(NA, n, n)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(sources)
colnames(p_mat) <- rownames(p_mat) <- colnames(sources)

for (i in 1:n) {
    for (j in 1:n) {
        test <- suppressWarnings(cor.test(sources[[i]], sources[[j]], method = "spearman"))
        cor_mat[i, j] <- test$estimate
        p_mat[i, j] <- test$p.value
    }
}

# png(file = "plots/microbial/sourcetracker/05_sources_correlation.png", width = 15, height = 10, units = "in", res = 600)
pheatmap::pheatmap(cor_mat,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         scale = "none",
         main = "Correlation between sources contributions")
# dev.off()

