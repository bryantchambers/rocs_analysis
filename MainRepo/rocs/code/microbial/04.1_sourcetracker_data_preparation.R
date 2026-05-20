library(gghighlight)
library(tidyverse)
library(biomformat)
library(phyloseq)
library(showtext)
library(ggh4x)
library(janitor)
library(ggvenn)
library(dendextend)

setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")
source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")



############################################
##### Create objects for sourcetracker #####
############################################

cdata <- read.table(file = "data/metadata_v5.tsv", sep = "\t", header = T) %>%
    filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))


# Load phyloseq object with filtered table
tax_data <- read_tsv("results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz") %>%
    filter(label %in% cdata$label) %>%
    filter(is_dmg == "Damaged")
    # filter(domain %in% c("d__Bacteria", "d__Archaea"))

tax_data_df <- tax_data %>%
    select(label, species, abundance) %>%
    group_by(label, species) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) |>
    as.data.frame() |>
    column_to_rownames("species")

taxonomy_df <- tax_data |>
    select(domain, phylum, class, order, family, genus, species) |>
    mutate(name = species) |>
    distinct() |>
    as.data.frame() |>
    column_to_rownames("name")

cdata_df <- cdata |>
    filter(label %in% colnames(tax_data_df)) |>
    mutate(name = label) |>
    distinct() |>
    as.data.frame() |>
    column_to_rownames("name") %>%
    mutate(SourceSink = "sink")

tax_data_phyloseq <- phyloseq::phyloseq(
    otu_table(tax_data_df, taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_df)),
    sample_data(cdata_df)
)


# Filter the phyloseq object
phyloseq_filt <- filt_physeq(tax_data_phyloseq, min_counts = 1000, ncounts = 3, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 1e-6)

cdata_sp <- sample_data(phyloseq_filt) 
st_sink <- phyloseq_filt


# # Source data
# source_cdata <- read_tsv("/projects/caeg/people/ngm902/apps/sourcetracker/source_data/med-biomes-download.txt")

# source_cdata_short <- source_cdata %>%
#     select(label = run_accession, biome, biome_class, biome_subclass) |>
#     # rename(Env = biome) |>
#     mutate(SourceSink = "source") |>
#     select(label, SourceSink, biome, biome_class, biome_subclass)

# source_data_lca <- read_tsv("/projects/caeg/people/ngm902/apps/sourcetracker/source_data/tp-lca.summary.tsv.gz") |>
#     filter(rank == "subspecies") |>
#     separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecie"), sep = ";") %>%
#     rename(nreads = n_reads)

# source_data_stats <- read_tsv("/projects/caeg/people/ngm902/apps/sourcetracker/source_data/tp-mapping-filtered.summary.tsv.gz") |>
#     filter(norm_entropy > 0.6, norm_gini < 0.4) |>
#     mutate(name = paste0("S__", reference)) |>
#     mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad))

# source_data <- source_data_lca |>
#     select(-nreads, -abundance) |>
#     inner_join(source_data_stats) |>
#     select(label, reference, abundance)

# # Glom at species level
# source_data_sp <- source_data_lca |> 
#     select(-nreads, -abundance) |>
#     inner_join(source_data_stats) %>%
#     # filter(domain %in% c("d__Bacteria", "d__Archaea")) |> 
#     group_by(label, species) |>
#     summarize(abundance = sum(abundance)) |>
#     ungroup()

# st_df <- source_data_sp |> 
#     pivot_wider(names_from = species, values_from = abundance, values_fill = 0) |>
#     as.data.frame() |>
#     column_to_rownames("label")

# st_taxonomy <- source_data_lca |>
#     select(domain, phylum, class, order, family, genus, species) %>%
#     distinct() %>%
#     inner_join(source_data_sp |> select(species) |> distinct()) |> 
#     select(domain, phylum, class, order, family, genus, species) |>
#     distinct() |>
#     mutate(reference = species) |>
#     distinct() |>
#     as.data.frame() |>
#     column_to_rownames("reference")

# st_cdata <- source_cdata_short |>
#     as.data.frame() |>
#     mutate(clabel = label) |>
#     column_to_rownames("clabel")

# st_source <- phyloseq(
#     otu_table(st_df, taxa_are_rows = FALSE),
#     tax_table(as.matrix(st_taxonomy)),
#     sample_data(st_cdata)
# )

# # Remove samples dominated by a single species
# to_rm <- microbiome::transform(st_source, "compositional") |>
#     speedyseq::psmelt() |>
#     as_tibble() |>
#     filter(Abundance > 0.9) |>
#     select(Sample) |>
#     distinct()

# st_source <- prune_samples(!(sample_names(st_source) %in% to_rm$Sample), st_source)

# st_source %>% saveRDS("/projects/caeg/people/ngm902/apps/repos/rocs/data/sourcetracker/phyloseq_sourcetracker_sources_wviruses.rds")



st_source <- readRDS("/projects/caeg/people/ngm902/apps/repos/rocs/data/sourcetracker/phyloseq_sourcetracker_sources_wviruses.rds")

cut_s <- st_sink |>
    sample_sums() |>
    min()

cut_s

st_source <- prune_samples(sample_sums(st_source) >= cut_s, st_source)
st_source <- prune_taxa(taxa_sums(st_source) > 0, st_source)



## How much of out sinks is on the sources?
st_sinks_sums <- st_sink %>%
    microbiome::transform("compositional") %>%
    prune_taxa(taxa_names(.) %in% taxa_names(st_source), .) %>%
    sample_sums() %>%
    enframe(name = "sample", value = "sum")

st_sinks_sums %>% 
    left_join(cdata, by = c("sample"="label")) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = sum))+
        geom_point()+
        geom_line()+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
        coord_flip()+
        facet_nested(~core, scales = "free", space = "free")


## Test removing samples in which shared taxa do not exceed a threshold
st_source_sums <- st_source %>%
    microbiome::transform("compositional") %>%
    prune_taxa(taxa_names(.) %in% taxa_names(st_sink), .) %>%
    sample_sums() %>%
    enframe(name = "sample", value = "sum")


# Opción A: aumentar número de barras usando `bins`
st_source_sums %>%
    ggplot(aes(x = sum)) +
        geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
        geom_histogram(bins = 100, color = "black", fill = "grey") +
        labs(x = "Relative abundance represented by the shared taxa", y = "Count")

st_source_sums %>%
    filter(sum >= 0.1) %>%
    nrow()

st_source <- st_source %>%
    prune_samples(sample_names(.) %in% st_source_sums[st_source_sums$sum >= 0.1,]$sample, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

st_source %>%
    sample_names() |>
    as_tibble() |>
    rename(label = value) |>
    head()

source_cdata <- read_tsv("/projects/caeg/people/ngm902/apps/repos/rocs/data/sourcetracker/med-biomes-download.txt") %>%
    mutate(biome_subclass = gsub(":|-| +", "_", biome_subclass)) %>%
    mutate(biome_subclass = tolower(biome_subclass)) %>%
    mutate(biome_subclass = gsub(" +", "_", biome_subclass))

source_cdata <- source_cdata %>%
    filter(run_accession %in% sample_names(st_source))
    

source_cdata %>% 
    count(biome) %>%
    arrange(desc(n)) %>%
    knitr::kable()

# |biome                                                        |  n|
# |:------------------------------------------------------------|--:|
# |root:Environmental:Aquatic:Marine                            | 73|
# |root:Environmental:Aquatic:Estuary                           | 49|
# |root:Environmental:Terrestrial:Soil                          | 31|
# |root:Environmental:Aquatic:Marine:Sediment                   | 25|
# |root:Environmental:Terrestrial:Soil:Agricultural             | 23|
# |root:Host-associated:Plants:Rhizosphere                      | 17|
# |root:Environmental:Aquatic:Marine:Oceanic                    | 13|
# |root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh | 10|
# |root:Environmental:Aquatic:Freshwater:Lotic:Sediment         |  6|
# |root:Environmental:Aquatic:Marine:Coastal                    |  6|
# |root:Environmental:Aquatic:Freshwater:Sediment               |  4|
# |root:Environmental:Terrestrial:Soil:Permafrost               |  3|
# |root:Environmental:Aquatic:Freshwater:Lake                   |  2|
# |root:Environmental:Terrestrial:Soil:Grasslands               |  2|
# |root:Host-associated:Plants                                  |  1|


# ## Standardize sample counts to the minimum
# st_ps <- merge_phyloseq(st_source, st_sink)

# st_ps |>
#     sample_sums() |>
#     summary()

# total <- min(sample_sums(st_ps))
# standf <- function(x, t = total) round(t * (x / sum(x)))
# st_ps <- transform_sample_counts(st_ps, standf)

# st_ps |>
#     sample_sums() |>
#     summary()
# st_biom <- biomformat::make_biom(
#     data = t(as((otu_table(st_ps, taxa_are_rows = FALSE)), "matrix")),
#     matrix_element_type = "int"
# )
# biomformat::write_biom(st_biom, biom_file = "results/microbial/sourcetracker/st_sp_min.biom")

# st_ps |>
#     speedyseq::psmelt() |>
#     as_tibble() |>
#     select(Sample, SourceSink, Env = biome_subclass) |>
#     distinct() |>
#     rename("#SampleID" = Sample) |>
#     write_tsv("results/microbial/sourcetracker/st_sp_min.map")


## Standardize sample counts to the median
st_ps <- merge_phyloseq(st_source, st_sink)

st_ps |>
    sample_sums() |>
    summary()

total <- median(sample_sums(st_ps))
standf <- function(x, t = total) round(t * (x / sum(x)))
st_ps <- transform_sample_counts(st_ps, standf)

st_ps |>
    sample_sums() |>
    summary()
st_biom <- biomformat::make_biom(
    data = t(as((otu_table(st_ps, taxa_are_rows = FALSE)), "matrix")),
    matrix_element_type = "int"
)
biomformat::write_biom(st_biom, biom_file = "results/microbial/sourcetracker/st_sp_median_wviruses.biom")

st_ps |>
    speedyseq::psmelt() |>
    as_tibble() |>
    select(Sample, SourceSink, Env = biome_subclass) |>
    distinct() |>
    rename("#SampleID" = Sample) |>
    write_tsv("results/microbial/sourcetracker/st_sp_median_wviruses.map")









############################################
##### Process sourcetracker outputs  #######
############################################

#### Output sourcetracker results
cdata <- read.table(file = "data/metadata_v5.tsv", sep = "\t", header = T)

source_cdata <- read_tsv("data/sourcetracker/med-biomes-download.txt") %>%
    mutate(biome_subclass = gsub(":|-| +", "_", biome_subclass)) %>%
    mutate(biome_subclass = tolower(biome_subclass)) %>%
    mutate(biome_subclass = gsub(" +", "_", biome_subclass))

header_line <- readLines("results/microbial/sourcetracker/st_sp_median_wviruses/mixing_proportions.txt", n = 1) %>% sub("^#", "", .)
stout <- read.table(file = "results/microbial/sourcetracker/st_sp_median_wviruses/mixing_proportions.txt", sep = "\t", header = FALSE, col.names = unlist(strsplit(header_line, "\t"))) %>% 
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
# colors[names(colors) == "Unknown"] <- "lightgrey"

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
            scale_y_sqrt(labels = scales::percent_format(), breaks = c(0, .05, .25, .50, .75))+
            labs(x = "", y = "Mean source proportion", fill = "Source")+
            # labs(x = "", y = "Mean source proportion", fill = "Source", title = "Mean source proportion by biome class and core")+
            theme(legend.position = "none")

png(file = "plots/microbial/sourcetracker/01_proportion_by_core.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()


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
        # labs(x = "", y = "Source proportion", fill = "Source", title = "Relative abundance of taxa assigned to each potential source") +
        scale_y_continuous(labels = scales::percent_format(), breaks = c(0, .2, .4, .6, .8)) +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4)) +
        facet_nested(~core, scales = "fixed", space = "fixed") +
        scale_fill_manual(values = colors) +
        coord_flip(ylim = c(0, 0.6))+
        guides(fill = guide_legend(reverse = TRUE))
        # theme(
            # legend.position = "bottom",
            # legend.direction = "horizontal",
            # legend.box = "horizontal",
            # aspect.ratio = 5/1)

png(file = "plots/microbial/sourcetracker/02_proportion_sample.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()



plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    filter(source %in% c("Polar", "Atlantic Ocean", "Mediterranean Sea")) %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "Age (y bp)", y = "Source proportion", fill = "Source")+
            scale_y_continuous(labels = scales::percent_format())+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
            facet_nested(~source*core, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()+
            guides(fill = guide_legend(nrow = 1)) +
            theme(
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.box = "horizontal",
                # aspect.ratio = 3/1
            )

# png(file = "results/plots/microbial/sourcetracker_03_mediterranean.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()


plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    filter(source == "Atlantic Ocean") %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "", y = "Source proportion", fill = "Source", title = "Relative abundance of taxa assigned to each potential source")+
            scale_y_continuous(labels = scales::percent_format())+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
            facet_nested(~core, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()+
            guides(fill = guide_legend(nrow = 1)) +
            theme(
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.box = "horizontal",
                # aspect.ratio = 3/1
            )
plotA


plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    filter(source == "Polar") %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "", y = "Source proportion", fill = "Source", title = "Relative abundance of taxa assigned to each potential source")+
            scale_y_continuous(labels = scales::percent_format())+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
            facet_nested(~core, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()+
            theme(
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.box = "horizontal",
                # aspect.ratio = 3/1
            )

# png(file = "results/plots/microbial/sourcetracker_03_polar.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()


plotA <- stout %>%
    left_join(cdata) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    filter(source %in% c("Brackish Water")) %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    filter(core == "ST13") %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "", y = "Source proportion", fill = "Source", title = "Relative abundance of taxa assigned to each potential source")+
            scale_y_continuous(labels = scales::percent_format())+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
            facet_nested(~core, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()
            # theme(
            #     legend.position = "bottom",
            #     legend.direction = "horizontal",
            #     legend.box = "horizontal",
            #     # aspect.ratio = 3/1
            # )

# png(file = "results/plots/microbial/sourcetracker_03_brackish.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()


plotA <- stout %>%
    left_join(cdata) %>%
    filter(grepl("Polar|Atlantic Ocean", source)) %>%
    select(label, y_bp, source, proportion) %>%
    # mutate(proportion = proportion + 1e-10) %>%
    pivot_wider(names_from = source, values_from = proportion, values_fill = 0) %>%
    clean_names() %>%
    mutate(total = atlantic_ocean + polar) %>%
    mutate(ratio = log2((atlantic_ocean + 1e-10) / (polar + 1e-10))) %>%
    # filter(total > 0.0) %>%
    pivot_longer(cols = c(atlantic_ocean, polar, ratio), names_to = "stat", values_to = "value") %>% 
    filter(!is.nan(value)) %>%
    filter(value != "Invalid Number") %>%
    left_join(cdata) %>%
    filter(stat == "ratio") %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = value)) +
        geom_area(stat = "identity", position = "stack", fill = "lightblue", alpha = 0.5) +
        geom_line(alpha = 0.3) +
        geom_hline(data = data.frame(stat = "ratio", y = 0), aes(yintercept = y), linetype = "dashed", color = "grey") +
        coord_flip() +
        facet_nested(~core, scales = "free", space = "fixed")+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
        labs(x = "Age (yp)", y = "Log2 ratio of Atlantic Ocean to Polar source proportion") +
        theme(aspect.ratio = 5/1)

png(file = "plots/microbial/sourcetracker/sourcetracker_04_atlantic_polar_ratio_s.png", width = 10, height = 10, units = "in", res = 600)
plotA
dev.off()






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
    geom_area(aes(y = ratio_pos), fill = "#FFC300") +
    geom_area(aes(y = ratio_neg), fill = "#154360") +
    geom_line(aes(y = ratio), color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    coord_flip() +
    facet_nested(~core, scales = "free", space = "fixed") +
    scale_x_reverse(
        labels = scales::label_number(scale_cut = scales::cut_short_scale()),
        breaks = seq(0, 5e06, 5e4)
    ) +
    labs(
        x = "Age (yBP)",
        y = "Log2 ratio of Atlantic Ocean to Polar source proportion"
    ) +
    # scale_fill_manual(values = c("#154360", "#FFC300"))+
    theme(aspect.ratio = 5 / 1)

png(file = "plots/microbial/sourcetracker/sourcetracker_04_atlantic_polar_ratio_s.png", width = 11, height = 10, units = "in", res = 600)
plotA
dev.off()

pdf(file = "plots/microbial/sourcetracker/sourcetracker_04_atlantic_polar_ratio_s.pdf", width = 11, height = 10)
plotA
dev.off()














plotA <- stout %>%
    left_join(cdata) %>%
    filter(grepl("Polar|Atlantic Ocean", source)) %>%
    select(label, y_bp, source, proportion) %>%
    # mutate(proportion = proportion + 1e-10) %>%
    pivot_wider(names_from = source, values_from = proportion, values_fill = 0) %>%
    clean_names() %>%
    mutate(total = atlantic_ocean + polar) %>%
    mutate(ratio = log2(atlantic_ocean/polar)) %>% 
    # filter(total > 0.0) %>%
    pivot_longer(cols = c(atlantic_ocean, polar, ratio), names_to = "stat", values_to = "value") %>%
    filter(!is.nan(value)) %>%
    left_join(cdata) %>%
    filter(stat == "ratio") %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = value)) +
        geom_point()+
        geom_line(alpha = 0.3) +
        geom_hline(data = data.frame(stat = "ratio", y = 0), aes(yintercept = y), linetype = "dashed", color = "grey") +
        coord_flip() +
        facet_nested(~core, scales = "free", space = "fixed")+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
        labs(x = "Age (yp)", y = "Log2 ratio of Atlantic Ocean to Polar source proportion") +
        theme(aspect.ratio = 5/1)

png(file = "plots/microbial/sourcetracker/sourcetracker_04_atlantic_polar_ratio.png", width = 11, height = 10, units = "in", res = 600)
plotA
dev.off()


ratios <- stout %>%
    left_join(cdata) %>%
    filter(grepl("Polar|Atlantic Ocean", source)) %>%
    select(label, source, proportion) %>%
    # mutate(proportion = proportion + 1e-10) %>%
    pivot_wider(names_from = source, values_from = proportion, values_fill = 0) %>%
    clean_names() %>%
    mutate(total = atlantic_ocean + polar) %>%
    mutate(ratio = log2(atlantic_ocean/polar)) %>% 
    # filter(total > 0.0) %>%
    pivot_longer(cols = c(atlantic_ocean, polar, ratio), names_to = "stat", values_to = "value") %>%
    filter(!is.nan(value)) %>%
    # left_join(cdata) %>%
    filter(stat == "ratio") %>%
    select(label, value)

plotA <- stout %>% 
    left_join(cdata) %>% 
    # filter(y_bp < 140000) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    # filter(biome_class %in% c("Environmental:Aquatic:Marine") | biome_class == "Unknown") %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    filter(core %in% c("GeoB25202_R2")) %>%
    # mutate(core = fct_relevel(core, "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "", y = "Source proportion", fill = "Source", title = "Relative abundance of taxa assigned to each potential source")+
            scale_y_continuous(labels = scales::percent_format())+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
            facet_nested(~core*source, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()

# png(file = "results/plots/microbial/sourcetracker_05_full_example_GeoB.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()

plotA <- stout %>% 
    left_join(cdata) %>% 
    # filter(y_bp < 140000) %>%
    group_by(source) %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    # filter(biome_class %in% c("Environmental:Aquatic:Marine") | biome_class == "Unknown") %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    mutate(source = fct_reorder(source, -mean_proportion)) %>%
    filter(core %in% c("ST13")) %>%
    # mutate(core = fct_relevel(core, "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source)) +
            geom_area(stat = "identity", position = "stack")+
            labs(x = "", y = "Source proportion", fill = "Source", title = "Relative abundance of taxa assigned to each potential source")+
            scale_y_continuous(labels = scales::percent_format())+
            scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
            facet_nested(~core*source, scales = "free", space = "fixed")+
            scale_fill_manual(values = colors)+
            coord_flip()

# png(file = "results/plots/microbial/sourcetracker_05_full_example_ST13.png", width = 15, height = 10, units = "in", res = 600)
plotA
# dev.off()

plotA <- cdata %>%
    filter(core == "GeoB25202_R2") %>%
    ggplot(., aes(x = y_bp, y = -mis)) +
        geom_point()+
        geom_line()+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
        coord_flip()+
        labs(x = "Age (y bp)", y = "d18O", title = "MIS (literature)")

plotB <- cdata %>%
    filter(core == "GeoB25202_R2") %>%
    ggplot(., aes(x = y_bp, y = temp)) +
        geom_point()+
        geom_line()+
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
        coord_flip()+
        labs(x = "", y = "Temperature (ºC)", title = "SST (in-house)")+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )

png(file = "results/plots/microbial/sst_mis.png", width = 8, height = 10, units = "in", res = 600)
egg::ggarrange(plotA, plotB, nrow = 1)
dev.off()


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


png(file = "results/plots/microbial/sourcetracker_05_sources_correlation.png", width = 15, height = 10, units = "in", res = 600)
# plot(2,4)
pheatmap::pheatmap(cor_mat,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         scale = "none",
         main = "Correlation between Sources")
dev.off()





#### Taxonomic profile for sourcetracker sources
tax_data <- read_tsv("results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz") %>%
    filter(label %in% cdata$label) %>%
    filter(is_dmg == "Damaged")


# Ruta a la carpeta que contiene los archivos
# folder <- "results/sourcetracker/st_sp_median/"
folder <- "results/microbial/sourcetracker/st_sp_median_wviruses/"

# Obtener lista de archivos .feature_table.txt
files <- list.files(folder, pattern = "\\.feature_table\\.txt$", full.names = TRUE)

# Crear lista vacía para guardar resultados
results_list <- list()
i <- 1
# Loop por cada archivo
for (i in seq_along(files)[1:20]) {
  file_path <- files[i]
  
  df <- read.table(file_path, sep = "\t", header = TRUE, check.names = FALSE) %>%
    rename("source" = 1) %>%
    column_to_rownames(var = "source") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    clean_names() %>%
    rownames_to_column("species") %>%
    pivot_longer(cols = -species, names_to = "source", values_to = "proportion") %>%
    mutate(
      source = gsub("^root_environmental_aquatic_", "", source),
      source = gsub("^root_environmental_terrestrial_soil_", "", source),
      source = gsub("^marine_", "", source),
      source = source %>%
        str_replace_all("_", " ") %>%
        str_to_lower() %>%
        str_to_title(),
      sample = basename(file_path) %>%
        str_remove("\\.feature_table\\.txt$")
    )
  
  # Guardar el resultado en la lista
  results_list[[i]] <- df
  print(paste0(i, "/", max(seq_along(files))))
}

# Unir todos los data.frames en uno solo
all_data <- bind_rows(results_list)
# write_tsv(all_data, "results/microbial/sourcetracker/st_sp_median/all_sources_long.tsv")


# all_data <- read_tsv("results/sourcetracker/st_sp_median/all_sources_long.tsv") %>%
all_data.b <- all_data %>%
    left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    group_by(species) %>%
    mutate(total_proportion = sum(proportion)) %>%
    ungroup() %>%
    filter(total_proportion > 0) %>%
    select(-total_proportion) %>%
    left_join(tax_data)


# SPECIES LEVEL
# 1. Crear matriz de abundancia (filas = source, columnas = phylum, valores = abundancia)
abund_matrix <- all_data.b %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus)) %>%
    group_by(source, domain, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    select(source, phylum, rel_proportion) %>%
    pivot_wider(names_from = phylum, values_from = rel_proportion, values_fill = 0) %>%
    column_to_rownames("source") %>%
    as.matrix()

# 2. Calcular distancia y hacer clustering
dist_mat <- vegan::vegdist(abund_matrix, method = "bray")
clust <- hclust(sqrt(dist_mat), method = "average")  # UPGMA

# 3. Obtener orden jerárquico de los sources
ordered_sources <- clust$labels[clust$order]

clust$labels <- gsub("North Pacific Subtropical Gyre", "North Pacific\nSubtropical Gyre", clust$labels)

dend <- as.dendrogram(clust)
# dend <- dendextend::rotate(dend)


size <- stout %>% 
    group_by(source) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
    filter(source %in% ordered_sources) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>%
    arrange(source) %>%
    mutate(source = gsub("North Pacific Subtropical Gyre", "North Pacific\nSubtropical Gyre", source))

# Función para añadir puntos al final de las ramas
add_points_to_leaves <- function(dend, size_df) {
    dendrapply(dend, function(node) {
        if(is.leaf(node)) {
            label <- attr(node, "label")
            idx <- which(size_df$source == label)
            if(length(idx) == 1) {
                attr(node, "nodePar") <- c(attr(node, "nodePar"),
                                           list(pch = 21,
                                                # cex = size_df$mean_proportion[idx] * 10,
                                                cex = -1/log10(size_df$mean_proportion[idx]),
                                                col = "black",
                                                bg = "#377EB8"))
            }
        }
        node
    })
}

dend <- add_points_to_leaves(dend, size)

# Ajusta los márgenes para que los textos del dendrograma no queden fuera
png(file = "results/plots/microbial/sourcetracker_05_1_species_contribution.png", width = 7, height = 15, units = "in", res = 600)
# op <- par(no.readonly = TRUE)
# par(mar = c(5, 4, 4, 10) + 0.1)  # aumenta el margen derecho
plot(dend, horiz = TRUE, leaflab = "none")  # remove tip labels
# plot(dend, horiz = TRUE)  # remove tip labels
# par(op)
dev.off()


plotB <- all_data.b %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(domain)) %>%
    # mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path")+
        theme(legend.position = "right")

png(file = "sourcetracker_05_2_genus_contribution.png", width = 15, height = 10, units = "in", res = 600)
plotB
dev.off()

plotB <- all_data %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(grepl("Sargasso", source)) %>% View()
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path")+
        theme(legend.position = "right")

png(file = "results/plots/microbial/sourcetracker_05_2_genus_contribution.png", width = 15, height = 10, units = "in", res = 600)
plotB
dev.off()


plotS <- all_data %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus, "; ", species)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(if_any(everything(), ~ grepl("Nitrosopumilaceae|JACPRH01", .))) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path", title = "Species level")+
        theme(legend.position = "none")

plotG <- all_data %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(family, "; ", genus)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(if_any(everything(), ~ grepl("Nitrosopumilaceae|JACPRH01", .))) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 2)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path", title = "Genus level")+
        theme(legend.position = "bottom")

png(file = "results/plots/microbial/sourcetracker_05_2_AOA.png", width = 20, height = 10, units = "in", res = 600)
egg::ggarrange(plotG, plotS, nrow = 1)
dev.off()


plotS <- all_data %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus, "; ", species)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(if_any(everything(), ~ grepl("Flavobacteriaceae", .))) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 2)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path", title = "Species level")+
        theme(legend.position = "none")

plotG <- all_data %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(family, "; ", genus)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(if_any(everything(), ~ grepl("Flavobacteriaceae", .))) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 3)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path", title = "Genus level")+
        theme(legend.position = "bottom")

png(file = "results/plots/microbial/sourcetracker_05_2_Flavobacteriaceae.png", width = 20, height = 10, units = "in", res = 600)
egg::ggarrange(plotG, plotS, nrow = 1)
dev.off()


plotS <- all_data %>% head
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("S") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus, "; ", species)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(if_any(everything(), ~ grepl("Pelagibacter", .))) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 2)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path", title = "Species level")+
        theme(legend.position = "none")

plotG <- all_data %>% 
    # filter(order == "o__Pelagibacterales") %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # left_join(tax_data) %>%
    mutate(phylum = paste0(family, "; ", genus)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.92)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    filter(if_any(everything(), ~ grepl("Flavobacteriaceae", .))) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 3)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path", title = "Genus level")+
        theme(legend.position = "bottom")

png(file = "results/plots/microbial/sourcetracker_05_2_Flavobacteriaceae.png", width = 20, height = 10, units = "in", res = 600)
egg::ggarrange(plotG, plotS, nrow = 1)
dev.off()

# plotB <- all_data %>% 
#     filter(source %in% c("Polar", "Atlantic Ocean")) %>%
#     mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus, "; ", species)) %>%
#     group_by(source, biome_class, phylum) %>%
#     summarise(proportion = sum(proportion), .groups = "drop") %>%
#     ungroup() %>%
#     # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
#     group_by(source) %>%
#     mutate(rel_proportion = proportion / sum(proportion)) %>%
#     ungroup() %>%
#     filter(source %in% unique(source[rel_proportion != 0])) %>%
#     group_by(phylum) %>%
#     mutate(total_abundance = sum(rel_proportion)) %>%
#     ungroup() %>%
#     mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
#     # filter(grepl("p__Thermoproteota", phylum)) %>% 
#     filter(total_abundance > quantile(total_abundance, 0.9)) %>%
#     # filter(total_abundance > . %>% pull(total_abundance) %>% unique() %>% sort() %>% slice(1:10)) %>% 
#     mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
#     ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
#         geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
#         coord_flip() +
#         # facet_nested(biome_class~., scales = "free", space = "free") +
#         scale_fill_manual(values = paired_genus) +
#         guides(fill = guide_legend(ncol = 1)) +
#         labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path")+
#         theme(legend.position = "right")

# #png(file = "results/plots/microbial/sourcetracker_05_2_species_contribution.png", width = 15, height = 10, units = "in", res = 600)
# plotB
# #dev.off()


plotB <- all_data %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    mutate(phylum = paste0(domain, "; ", phylum)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    # mutate(rel_proportion = proportion) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    filter(total_abundance > quantile(total_abundance, 0.5)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(.~phylum, scales = "free", space = "fixed") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "", y = "Family relative contribution", fill = "Taxonomic path")+
        theme(legend.position = "right")

png(file = "results/plots/microbial/sourcetracker_05_2_phylum_contribution.png", width = 15, height = 10, units = "in", res = 600)
plotB
dev.off()

plotB <- all_data %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    filter(total_abundance > quantile(total_abundance, 0.85)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "", y = "Family relative contribution", fill = "Taxonomic path")+
        theme(legend.position = "right")

png(file = "results/plots/microbial/sourcetracker_05_2_family_contribution.png", width = 15, height = 10, units = "in", res = 600)
plotB
dev.off()





## Taxonomic association with sources
unknown_association <- all_data %>% 
    # filter(source %in% c("Polar", "Unknown")) %>% 
    filter(biome_class %in% c("Environmental:Aquatic:Marine") | source == "Unknown") %>%
    # filter(source %in% c("Polar", "Atlantic Ocean", "Mediterranean Sea", "Unknown")) %>% 
    mutate(source_edited = ifelse(source == "Unknown", "Unknown", "Other")) %>%
    group_by(species, source_edited, sample, biome_class) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%

    group_by(species, sample) %>%
    mutate(sum_proportion = sum(proportion)) %>%
    ungroup() %>%
    filter(sum_proportion > 0) %>%
    select(-sum_proportion) %>% 
    
    group_by(species, sample) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>% 
    
    group_by(species, source_edited) %>%
    summarise(mean_rel_proportion = mean(rel_proportion), median_rel_proportion = median(rel_proportion), .groups = "drop")  %>%

    select(species, source_edited, median_rel_proportion) %>%
    pivot_wider(names_from = source_edited, values_from = median_rel_proportion, values_fill = 0) %>%
    janitor::clean_names() 

unknown_association %>%
    pivot_longer(cols = c(unknown, other), names_to = "source", values_to = "median_rel_proportion") %>%
    ggplot(aes(x = source, y = median_rel_proportion)) +
        geom_boxplot()+
        geom_violin(fill = "lightgrey", alpha = 0.5, width = 0.5)+
        coord_flip()

unknown_taxa <- unknown_association %>% 
    filter(other >= 0.5) %>%
    pull(species)



## Venn diagram polar and atlantic ocean + mediterranean sea
# Crear una lista de especies por source
df_species_by_source <- all_data %>%
    filter(species %in% unknown_taxa) %>%
    # filter(source %in% c("Polar", "Unknown")) %>% 
    filter(source %in% c("Polar", "Atlantic Ocean", "Mediterranean Sea")) %>% 
    # filter(source %in% c("Polar", "Atlantic Ocean")) %>% 
    mutate(source_edited = ifelse(source == "Polar", "Polar", "Atlantic Ocean")) %>%
    group_by(species, source_edited, sample, biome_class) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    # filter(source %in% as.character(stout %>% filter(biome_class == "Environmental:Aquatic:Marine") %>% pull(source) %>% unique())) %>%
    filter(proportion > 0) %>%
    # View()
    group_by(source_edited) %>%
    summarise(species = list(unique(species))) %>%
    deframe()

# Graficar diagrama de Venn (máx 4-5 sources para visualización clara)
selected_sources <- names(df_species_by_source)[1:length(df_species_by_source)]
venn_list <- df_species_by_source[selected_sources]

## Association strength between Polar and Atlantic Ocean
# all_data <- all_data %>% 
#     left_join(stout %>% select(source, biome_class) %>% distinct()) 

tax_association <- all_data %>%
    filter(species %in% unknown_taxa) %>%
    # filter(source %in% c("Polar", "Unknown")) %>% 
    filter(source %in% c("Polar", "Atlantic Ocean", "Mediterranean Sea")) %>% 
    # filter(source %in% c("Polar", "Atlantic Ocean")) %>% 
    group_by(species) %>%
    mutate(total_proportion = sum(proportion)) %>%
    ungroup() %>%
    filter(total_proportion > 0) %>%
    select(-total_proportion) %>%
    mutate(source_edited = ifelse(source == "Polar", "Polar", "Atlantic Ocean")) %>%
    group_by(species, source_edited, sample, biome_class) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    
    group_by(species, sample) %>%
    mutate(sum_proportion = sum(proportion)) %>%
    ungroup() %>%
    filter(sum_proportion > 0) %>%
    select(-sum_proportion) %>% 
    
    group_by(species, sample) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>% 
    
    group_by(species, source_edited) %>%
    summarise(mean_rel_proportion = mean(rel_proportion), median_rel_proportion = median(rel_proportion), .groups = "drop")  %>%

    select(species, source_edited, median_rel_proportion) %>%
    pivot_wider(names_from = source_edited, values_from = median_rel_proportion, values_fill = 0) %>%
    janitor::clean_names() %>%
    mutate(ratio = log2(atlantic_ocean/polar))
    

## Plotting
library(patchwork)

p1 <- ggvenn(venn_list, fill_color = c("#E41A1C", "#377EB8"), stroke_size = 0.5, set_name_size = 4)

p2 <- tax_association %>%
    # filter(polar > 0) %>%
    ggplot(aes(x = polar)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 50, color = "black", alpha = 0.7) +
    scale_fill_gradient2(low = "#E41A1C", mid = "white", high = "#377EB8", midpoint = 0.5) +
    labs(x = "Association to polar source", y = "Frequency", fill = "", title = "") +
    theme(legend.position = "none", aspect.ratio = 1/3)

p1 / p2






dmg_ssp <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp_selected_ss.tsv.gz") %>%
    filter(selected == "yes") #%>%
    # group_by(label) %>%
    # mutate(abundance = abundance / sum(abundance)) %>%
    # ungroup()

# dmg_ssp <- dmg_ssp %>%
#     # filter(selected == "yes", sediment_specific == "non-specific") %>%
#     group_by(label) %>%
#     mutate(abundance = abundance / sum(abundance)) %>%
#     ungroup()

dmg_ssp %>%
    mutate(type = case_when(
        species %in% as.character(tax_association %>% filter(polar >= 0.5) %>% pull(species) %>% unique()) ~ "polar",
        species %in% as.character(tax_association %>% filter(atlantic_ocean >= 0.5) %>% pull(species) %>% unique()) ~ "atlantic",
        TRUE ~ "other"
    )) %>%
    group_by(type, label) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    complete(type, label, fill = list(abundance = 0)) %>%
    left_join(cdata) %>%
    filter(type != "other") %>%
    mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = abundance, fill = type)) +
        geom_area(stat = "identity", position = "stack") +
        labs(x = "", y = "Abundance", color = "Phylum", title = "Abundance of taxa associated with polar source") +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4)) +
        facet_nested(.~core*type, scales = "free") +
        scale_fill_manual(values = paired_genus) +
        # coord_flip(xlim = c(140000, 0))
        coord_flip()

dmg_ssp %>%
    filter(core == "ST13") %>%
    mutate(type = case_when(
        species %in% as.character(tax_association %>% filter(polar >= 0.75) %>% pull(species) %>% unique()) ~ "polar",
        species %in% as.character(tax_association %>% filter(atlantic_ocean >= 0.75) %>% pull(species) %>% unique()) ~ "atlantic",
        TRUE ~ "other"
    )) %>%
    group_by(phylum, type, label) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    complete(phylum, type, label, fill = list(abundance = 0)) %>%
    left_join(cdata) %>%
    # filter(type != "other") %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    # filter(total_abundance > quantile(total_abundance, 0.5)) %>%
    # mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    mutate(type = fct_relevel(type, "atlantic", "polar", "other")) %>%
    ggplot(., aes(x = y_bp, y = abundance, fill = phylum)) +
        geom_area(stat = "identity", position = "stack") +
        labs(x = "", y = "Abundance", color = "Phylum", title = "Abundance of taxa associated with polar source") +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4)) +
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        facet_nested(.~core*type, scales = "free") +
        guides(fill = guide_legend(ncol = 1))+
        scale_fill_manual(values = paired_genus) +
        coord_flip(xlim = c(140000, 0))


dmg_ssp %>% # pull(phylum) %>% unique()
    group_by(label) %>%
    mutate(abundance = abundance / sum(abundance)) %>%
    ungroup() %>%
    filter(core == "ST13") %>%
    # filter(phylum == "p__Verrucomicrobiota") %>% 
    filter(order == "o__Rhodobacterales") %>%
    # filter(family == "f__Pseudoalteromonadaceae") %>%
    mutate(type = case_when(
        species %in% as.character(tax_association %>% filter(polar >= 0.75) %>% pull(species) %>% unique()) ~ "polar",
        species %in% as.character(tax_association %>% filter(atlantic_ocean >= 0.75) %>% pull(species) %>% unique()) ~ "atlantic",
        TRUE ~ "other"
    )) %>%
    mutate(taxa = paste(order, family, genus, sep = "; ")) %>%
    group_by(taxa, type, label) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    complete(taxa, type, label, fill = list(abundance = 0)) %>%
    left_join(cdata) %>%
    # filter(type != "other") %>%
    # filter(type != "polar") %>%
    group_by(taxa) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    filter(total_abundance > 0) %>%
    mutate(taxa = fct_reorder(taxa, -total_abundance)) %>%
    # filter(total_abundance > quantile(total_abundance, 0.5)) %>%
    # mutate(core = fct_relevel(core, "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    ggplot(., aes(x = y_bp, y = abundance, fill = taxa)) +
        geom_area(stat = "identity", position = "stack") +
        labs(x = "", y = "Abundance", color = "Phylum", title = "Abundance of taxa associated with polar source") +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4)) +
        # scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        facet_nested(.~core*type, scales = "free") +
        guides(fill = guide_legend(ncol = 1))+
        scale_fill_manual(values = paired_genus) +
        coord_flip(xlim = c(140000, 0))

dmg_ssp %>%
    select(domain, phylum, class, order, family, genus, species) %>%
    distinct() %>%
    arrange(domain, phylum, class, order, family, genus, species) %>%
    filter(species %in% as.character(tax_association %>% filter(atlantic_ocean >= 0.75) %>% pull(species) %>% unique())) %>%
    print(n = 100)

dmg_ssp %>%
    mutate(type = case_when(
        species %in% as.character(tax_association %>% filter(polar >= 0.75) %>% pull(species)) ~ "polar",
        species %in% as.character(tax_association %>% filter(atlantic_ocean >= 0.75) %>% pull(species)) ~ "atlantic",
        TRUE ~ "other"
    )) %>% 
    group_by(type, phylum, label) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    complete(type, phylum, label, fill = list(abundance = 0)) %>%
    left_join(cdata) %>% 
    filter(type != "other") %>%
    # filter(type == "atlantic") %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.95)) %>%
    # left_join(tax_data) %>%
    filter(phylum == c("p__SAR324", "p__Cyanobacteria")) %>%
    mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
    # filter(core == "GeoB25202_R1") %>%
    filter(abundance != 0)  %>%
    ggplot(., aes(x = y_bp, y = abundance, color = phylum)) +
        # geom_area(stat = "identity", position = "stack") +
        geom_point()+
        # geom_line()+
        labs(x = "", y = "Abundance", color = "Phylum", title = "Abundance of taxa associated with polar source") +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4)) +
        facet_nested(.~core, scales = "free") +
        guides(fill = guide_legend(ncol = 1)) +
        scale_fill_manual(values = paired_genus) +
        coord_flip()
 

dmg_ssp %>% 
    filter(if_any(c(domain, phylum, class, order, family, genus, species, subspecies), ~ grepl("g__Mazuvirus", .))) %>%
    group_by(genus, y_bp, core) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    ggplot(., aes(x = y_bp, y = abundance, color = genus)) +
        # geom_area(stat = "identity", position = "stack") +
        geom_point()+
        # geom_line()+
        labs(x = "", y = "Abundance", color = "Phylum", title = "Abundance of taxa associated with polar source") +
        scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4)) +
        facet_nested(.~core, scales = "free") +
        guides(fill = guide_legend(ncol = 1)) +
        scale_fill_manual(values = paired_genus) +
        coord_flip()

dmg_ssp %>%
    select(domain, phylum, class, order, family, genus, species, subspecies) %>%
    distinct() %>%
    arrange(domain, phylum, class, order, family, genus, species, subspecies) %>%
    filter(domain == "d__Viruses") %>%
    pull(genus) %>% unique()
    View()






dmg_ssp <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp_selected_ss.tsv.gz") %>%
    filter(selected == "yes")


data_wide <- dmg_ssp %>%
    filter(core == "GeoB25202_R1") %>%
    group_by(label) %>%
    mutate(abundance = abundance / sum(abundance)) %>%
    ungroup() %>%
    group_by(species) %>%
    mutate(count = n_distinct(label)) %>%
    ungroup() %>%
    filter(count > 10) %>%
    group_by(label, species) %>% 
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
    left_join(ratios %>% rename(valueeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee = value)) %>%
    column_to_rownames("label")

dim(data_wide)

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

# --- Corrección por múltiples comparaciones (FDR / Benjamini-Hochberg) ---
# Convertimos p.mat a vector, corregimos, y lo volvemos a matriz
pvals_vector <- p.mat[upper.tri(p.mat)]
padj_vector <- p.adjust(pvals_vector, method = "BH")

p.mat.adj <- matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
p.mat.adj[upper.tri(p.mat.adj)] <- padj_vector
p.mat.adj[lower.tri(p.mat.adj)] <- t(p.mat.adj)[lower.tri(p.mat.adj)]
diag(p.mat.adj) <- 0
colnames(p.mat.adj) <- colnames(p.mat)
rownames(p.mat.adj) <- rownames(p.mat)

dist_mat <- as.dist(1 - M)

hc <- hclust(dist_mat, method = "ward.D2")

old_par <- par(no.readonly = TRUE) # Guarda la configuración actual
par(mar = c(5, 20, 4, 4)) # Aumenta el margen derecho (último valor)
# plot(as.dendrogram(hc), horiz = TRUE, main = "Taxa clustering by correlation", xlab = "", sub = "")
dendextend::plot_horiz.dendrogram(hc, side = TRUE)
par(old_par) # Restaura la configuración original


# --- Heatmap reordenado según clustering ---
M_reordered <- M[rev(hc$order), rev(hc$order)]
p.mat_reordered <- p.mat.adj[rev(hc$order), rev(hc$order)]

corrplot::corrplot(M_reordered, method = "color", order = "original", 
  type = "upper", p.mat = p.mat_reordered, sig.level = 0.01, insig = "blank",
  tl.srt = 60, tl.cex = 0.7, tl.col = "black") # Reduce el tamaño de letra a 80% y texto negro

















maintaxa <- read.csv("/projects/caeg/people/ngm902/rocs/v3/main_references.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  select(-tax) %>%
  left_join(dmg_ssp %>% select(phylum, domain, class, order, family, genus, species, subspecies) %>% distinct())

dmg_ssp %>%
    mutate(type = case_when(
        species %in% as.character(tax_association %>% filter(polar >= 0.75) %>% pull(species)) ~ "polar",
        species %in% as.character(tax_association %>% filter(atlantic_ocean >= 0.75) %>% pull(species)) ~ "atlantic",
        TRUE ~ "other"
    )) %>%
    filter(type != "other") %>%
    group_by(phylum, class, type) %>%
    summarise(abundance = mean(abundance), .groups = "drop") %>%
    ggplot(., aes(x = type, y = abundance, fill = paste(phylum, class, sep = "; "))) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "", y = "Abundance", color = "Phylum", title = "Abundance of taxa associated with polar source") +
        scale_fill_manual(values = paired_genus) +
        theme(legend.position = "right") +
        coord_flip() +
        guides(fill = guide_legend(ncol = 1))
    
dmg_ssp %>%
    filter(subspecies == "S__GCA_905612295.1") %>%
    ggplot(., aes(x = y_bp, y = abundance, fill = paste(domain, phylum, class, order, family, genus, subspecies, sep = "; ")) ) +
		geom_area(stat = "identity", position = "stack") +
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
		coord_flip() +
		facet_nested(~core, scales = "free")+
		scale_fill_manual(values = paired_genus) +
		guides(fill = guide_legend(ncol = 1, title = "Taxa")) +
		labs(fill = "Taxa")




all_data %>% 
    left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    # filter(biome_class %in% c("Environmental:Aquatic:Marine", "Unknown")) %>%
    # filter(biome_class %in% c("Environmental:Aquatic:Marine", "Unknown")) %>%
    filter(source %in% c("Polar", "Atlantic Ocean")) %>%
    left_join(tax_data) %>%
    mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus, "; ", species)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>%
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    filter(source %in% unique(source[rel_proportion != 0])) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    # filter(grepl("p__Thermoproteota", phylum)) %>% 
    filter(total_abundance > quantile(total_abundance, 0.995)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path")+
        theme(legend.position = "right")





spec_index <- all_data %>%
    filter(source %in% as.character(stout %>% filter(biome_class == "Environmental:Aquatic:Marine") %>% pull(source) %>% unique())) %>%
    # filter(source %in% as.character(stout %>% filter(biome_class != "Unknown") %>% pull(source) %>% unique())) %>%
    group_by(, sample) %>%
    mutate(total_sp = sum(proportion, na.rm = TRUE)) %>%
    filter(total_sp > 0) %>%  
    mutate(rel_weight = proportion / total_sp) %>%
    ungroup() %>%
    group_by(species, source) %>%
    summarise(
    spec_index = mean(rel_weight, na.rm = TRUE),
    n_presences = n(),
    .groups = "drop"
    ) %>%
    arrange(desc(spec_index)) %>%
    left_join(tax_data)

sediment_specific_taxa <- spec_index %>%
    filter(source == "Sediment" & spec_index >= 0.95) %>% pull(phylum) %>% unique()
    pull(species) %>%
    append(c("s__DG-33", "s__DG-33 sp001515185", "s__DG-33 sp004375695", "s__DG-33 sp004375695")) # Hadarchaeota species that are sediment specific (references: S__3300024429_35, S__3300024433_66, S__3300001753_19, S__GCA_004375695.1, S__GCA_009619155.1)



dmg_ssp %>%
    group_by(domain, phylum, class, order, family, genus, species, subspecies) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    arrange(desc(abundance)) %>%
    filter(species == "s__Nitrosopumilus sp013203245") %>%
    print(n = 50)




p__Aerophobota
p__UBA6262
p__Hadarchaeota
p__Atribacterota
p__Asgardarchaeota
p__UBP18
p__KSB1
g__DRGT01 # Nitrosopulimaceae



dmg_ssp <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp_selected.tsv.gz")
dmg_ssp %>%    
    mutate(sediment_specific = case_when(
        phylum %in% c("p__Aerophobota", "p__UBA6262", "p__Hadarchaeota", "p__Atribacterota", "p__Asgardarchaeota", "p__UBP18", "p__KSB1") ~ "specific",
        genus %in% c("g__DRGT01") ~ "specific",
        # species %in% c("s__Nitrosopumilus sp013203245") ~ "specific",
        TRUE ~ "non-specific")) %>%

    # mutate(sediment_specific = ifelse(species %in% sediment_specific_taxa, "specific", "non-specific")) %>%
    write_tsv("results/microbial/taxonomy/dmg-summary-ssp_selected_ss.tsv.gz")

    






spec_index <- all_data %>%
    # filter(source != "Unknown") %>%
    mutate(source_type = ifelse(source %in% as.character(stout %>% filter(biome_class == "Environmental:Aquatic:Marine" & source != "Unknown") %>% pull(source) %>% unique()), "Marine", "Other")) %>%
    
    group_by(species, sample, source_type) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    
    group_by(species, sample) %>%
    filter(any(proportion != 0)) %>%
    ungroup() %>%

    group_by(species, sample) %>%
    mutate(relative_proportion = proportion/ sum(proportion)) %>%
    ungroup() %>% 
    
    group_by(species, source_type) %>%
    summarise(n = n_distinct(sample), mean_proportion = mean(relative_proportion), median_proportion = median(relative_proportion), .groups = "drop") %>%
    filter(source_type == "Marine") %>%
    rename(n_samples_st = n, mean_proportion_marine = mean_proportion, median_proportion_marine = median_proportion) %>%
    select(-source_type)







# Atlantic ocean vs Polar
# Atlantic ocean:
# Higher contribution of Talassarchaeaceae (genera MGIIB-02 and MGIIB-01) and Pelagicabter, 
# also present in Polar but at lower levels.
# High contributio of Gammaproteobacter (order SAR86; family D2472; genus D2472).
#
# Polar:
# Higher contribution of Nitrosopumilus compared to Atlantic Ocean.
# Arctic96AD-7 (order SAR324; family NAC60-12), polaribacter (family flavobacteriaceae)
# and other flavobacteriales abundant in polar but not in Atlantic Ocean.







dmg_ssp <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp_selected_ss.tsv.gz")
tax_data <- dmg_ssp %>%
    select(domain, phylum, class, order, family, genus, species, subspecies) %>%
    distinct()

dmg_ssp %>%
    group_by(species) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    left_join(spec_index) %>%
    arrange(desc(abundance)) %>%
    filter(mean_proportion_marine != "Invalid Number") %>%
    ggplot(., aes(x = abundance, y = mean_proportion_marine)) +
        geom_point()+
        geom_hline(yintercept = 0.0005, linetype = "dashed", color = "red")+
        scale_x_sqrt()+
        scale_y_log10()


dmg_ssp %>%
    filter(selected == "yes") %>%
    group_by(subspecies) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    left_join(tax_data) %>%
    left_join(spec_index) %>%
    filter(mean_proportion_marine != "Invalid Number") %>%
    arrange(desc(mean_proportion_marine)) %>%
    View()
    
marine <- spec_index %>%
    filter(mean_proportion_marine > 0.0005 | mean_proportion_marine == "Invalid Number")

dmg_ssp %>%
    filter(selected == "yes") %>%
    filter(core == "GeoB25202_R1") %>%
    mutate(marine_specific = ifelse(species %in% marine$species, "marine", "non-marine")) %>%
    group_by(phylum, marine_specific, label) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    complete(phylum, marine_specific, label, fill = list(abundance = 0)) %>%
    left_join(cdata) %>%
    group_by(phylum) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    # filter(grepl("Aerop", phylum)) %>%
    # mutate(phylum = fct_reorder(phylum, -total_abundance)) %>% 
    ggplot(., aes(x = y_bp, y = abundance, fill = phylum)) +
        geom_area(stat = "identity", position = "stack")+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
		coord_flip()+
		facet_nested(~core*marine_specific, scales = "fixed")+
		scale_fill_manual(values = paired_genus) +
		guides(fill = guide_legend(ncol = 1, title = "Species", reverse = TRUE))+
		theme(legend.position = "right")+
        labs(title = "..")
		# labs(title = paste0("Core ", target_core), x = "Age (y bp)", y = "Abundance") +

dmg_ssp %>%
    filter(selected == "yes") %>%
    # mutate(marine_specific = ifelse(species %in% marine$species, "marine", "non-marine")) %>%
    group_by(phylum, class, order, family, genus, species, subspecies) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    left_join(spec_index) %>%
    mutate(is_marine = ifelse(mean_proportion_marine > 0.0005 | is.na(mean_proportion_marine), "marine", "non-marine")) %>%
    


tax_data %>%
    filter(species == "s__SOIV01 sp011049545")

all_data %>%
    filter(grepl("JAEHKU01", species)) %>%
    filter(proportion > 0) %>%
    group_by(species, sample) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    # group_by(species, source) %>%
    # summarise(rel_proportion = mean(rel_proportion)) %>%
    View()

tax_data %>%
    filter(grepl("Nitrosopumilaceae", family)) %>%
    print(40)



all_data %>%
    left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    filter(biome_class %in% c("Environmental:Aquatic:Marine", "Unknown")) %>%
    # filter(biome_class %in% c("Environmental:Aquatic:Marine")) %>%
    left_join(tax_data) %>%
    filter(grepl("Nitrosopumilaceae", family)) %>% 
    # mutate(phylum = paste0(domain, "; ", phylum, "; ", class, "; ", order, "; ", family, "; ", genus, "; ", species)) %>%
    mutate(phylum = paste0(genus, "; ", species)) %>%
    group_by(source, biome_class, phylum) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    ungroup() %>% 
    # left_join(stout %>% select(source, biome_class) %>% distinct()) %>%
    group_by(source) %>%
    mutate(rel_proportion = proportion / sum(proportion)) %>%
    ungroup() %>% 
    filter(source %in% unique(source[rel_proportion != 0])) %>% 
    group_by(phylum) %>%
    mutate(total_abundance = sum(rel_proportion)) %>%
    ungroup() %>%
    # filter(total_abundance > quantile(total_abundance, 0.95)) %>%
    mutate(phylum = fct_reorder(phylum, -total_abundance)) %>%
    mutate(source = factor(source, levels = ordered_sources)) %>% # 4. Reordenar niveles del factor en tu df original
    # filter(source %in% c("Atlantic Ocean", "Polar", NA, "Sediment")) %>%
    # filter(source %in% c(NA, "Sediment")) %>%
    ggplot(., aes(x = source, y = rel_proportion, fill = phylum)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
        coord_flip() +
        # facet_nested(biome_class~., scales = "free", space = "free") +
        scale_fill_manual(values = paired_genus) +
        guides(fill = guide_legend(ncol = 2)) +
        labs(x = "", y = "Genera relative contribution", fill = "Taxonomic path")+
        theme(legend.position = "right")


dmg_ssp %>%
    group_by(species) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    arrange(desc(abundance)) %>%
    head(20)
    filter(grepl("sp013203245", species)) %>%
dmg_ssp %>%
    head()


