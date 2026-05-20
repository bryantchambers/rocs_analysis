library(tidyverse)
library(knitr)
library(readxl)
library(measurements)
library(ggmap)
library(ggrepel)
library(lvplot)
library(ggh4x)
library(phyloseq)

setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")

source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")
source("/projects/caeg/people/ngm902/scripts/dmg.R")
source("/projects/caeg/people/ngm902/scripts/get-metadata.R")

# httpgd::hgd()
# httpgd::hgd_browse()

# Load labels data
cdata <- read_tsv(file = "data/metadata_v5.tsv")
cdata_control <- cdata %>% filter(core %in% c("ExrNTC", "ExrPTC", "LibNTC", "LibPTC", "NEGATIVE")) %>% select(label, core, flowcell)


##########################################################################
####### Control samples plotting and references removal ##################
##########################################################################


tax_data <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")

tax_data_control <- tax_data %>%
	filter(core %in% c("ExrNTC", "ExrPTC", "LibNTC", "LibPTC", "NEGATIVE"))

# tax_data_control %>%
# 	ggplot(., aes(x = A_b, y = n_reads_ref, color = paste(domain, phylum)))+
# 		geom_point()+
# 		facet_wrap(~label)+
# 		scale_color_manual(values = paired_genus)

plotA <- tax_data_control %>%
	select(-core, -flowcell) %>%
	mutate(tax = paste(domain, phylum, class, order)) %>%
	group_by(tax) %>%
	mutate(reads_total = sum(n_reads_ref)) %>%
	ungroup() %>%
	mutate(tax = ifelse(reads_total < 1000, "Other (total nº reads < 1000)", tax)) %>%
	mutate(tax = fct_reorder(tax, -reads_total)) %>%
	left_join(cdata %>% select(label, core, flowcell)) %>%
	ggplot(., aes(x = label, y = A_b))+
		geom_point(aes(color = tax, size = n_reads_ref))+
		scale_color_manual(values = paired_genus)+
		labs(x = "Library", y = "Damage (A_b)", color = "Taxonomic path (domain, phylum, class, order)", size = "Number of reads")+
		guides(color = guide_legend(ncol = 1))+
		facet_nested(core~., scales = "free_y", space = "free_y")+
		coord_flip()

png(file = "plots/microbial/preliminary/control_libraries.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()

control_references <- tax_data_control %>%
	filter(!core %in% c("ExrPTC", "LibPTC") & label != "LV7008867195") %>%
	group_by(domain, phylum, class, order, family, genus, species, subspecies) %>%
	summarise(n = n(), n_reads_ref = sum(n_reads_ref), .groups = "drop") %>%
	arrange(desc(n_reads_ref))

# Let's filter out the control-specific references, and species which appear in more than 10 control samples 
tax_data_dmg_filtered <- tax_data_dmg %>%
	filter(
		!subspecies %in% control_references$subspecies & 
		!species %in% control_references[control_references$n >= 10,]$species &
		species != "s__Rhodococcus erythropolis")

# Remove low-complexity samples
label_stats <- tax_data_dmg_filtered %>%
	filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
	group_by(label) %>%
	summarise(n = n(), n_reads = sum(n_reads), .groups = "drop")

tax_data_dmg_filtered <- tax_data_dmg_filtered %>% 
	filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>% 
	filter(label %in% label_stats[label_stats$n_reads >= 1e4 & label_stats$n >= 30,]$label) 

tax_data_dmg_filtered %>%
	write_tsv("results/microbial/taxonomy/dmg-summary-ssp_clean.tsv.gz")