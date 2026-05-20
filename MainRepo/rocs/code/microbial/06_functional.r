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




# Load cdata
cdata <- read.table(file = "data/metadata_v5.tsv", sep = "\t", header = T) %>%
    filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))

# Load filtered table
tax_data <- read_tsv("results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz") %>%
    filter(label %in% cdata$label) %>%
    # filter(is_dmg == "Damaged") %>%
    filter(domain %in% c("d__Bacteria", "d__Archaea"))

# Load kegg data
kegg_data <- read_tsv("data/functional/kegg-modules-summary-rocs.tsv.gz", col_types = cols(gene_caller_ids_in_module="c"))

# Load sap data
sap_data <- read_tsv("data/functional/sap-predictions-rocs.tsv.gz")






# Relative abundance excluding viruses
table <- tax_data |> 
    filter(domain %in% c("d__Bacteria", "d__Archaea")) |> 
    mutate(reference = gsub("S__", "", subspecies)) %>%
    mutate(across(c(domain, phylum, class, order, family, genus, species), ~ str_remove(., "^[dpcfogs]__")))


cytochrome <- kegg_data %>% 
    filter(module == "M00155" & pathwise_module_completeness >= 1) %>%
    pull(reference)


