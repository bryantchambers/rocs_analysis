# Set up
library(tidyverse)
library(knitr)
library(readxl)
library(measurements)
library(ggh4x)
library(ggrepel)
library(dplyr)

source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")
setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")

cdata <- read_tsv("data/metadata_v5.tsv")

cdata %>%
	select(label, core, depth_in_core_cm) %>%
	mutate(depth = 
		case_when(
			core == "ST5" ~ -1911,
			core == "ST8" ~ -1908,
			core == "ST13" ~ -1731,
			core %in% c("GeoB25202_R1", "GeoB25202_R2") ~ NA,
			TRUE ~ NA_real_
		)) %>%
	mutate(data_collected = case_when(
		core == "ST5" ~ "03.03.22", 
		core == "ST8" ~ "09.09.21", 
		core == "ST13" ~ "08.12.21", 
		core %in% c("GeoB25202_R1", "GeoB25202_R2") ~ NA,
		TRUE ~ NA_character_
	)) %>%
	mutate(lat = case_when(
		core == "ST5" ~ "61_56_88", 
		core == "ST8" ~ "64_33_10",
		core == "ST13" ~ "63_15_02",
		core %in% c("GeoB25202_R1", "GeoB25202_R2") ~ NA,
		TRUE ~ NA_real_
	)) %>%
	mutate(lon = case_when(
		core == "ST5" ~ "19'29.70''", 
		core == "ST8" ~ "29'32.522",
		core == "ST13" ~ "28'15.24",
		core %in% c("GeoB25202_R1", "GeoB25202_R2") ~ NA,
		TRUE ~ NA_real_
	)) %>% head()

lat = c(, , , 63.5),
lon = c(-19.5, -29.5, -28.3, -28.5)

###################################################################
######## Sea Surfate Temeperature #################################
###################################################################

plotA1 <- cdata %>%
	filter(core != "GeoB25202_R2") %>%	
	mutate(core = gsub("_R1", "", core)) %>%
	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202"))) %>%
	filter(!is.na(temp_method)) %>%
	ggplot(., aes(x = y_bp, y = temp, color = temp_method))+
		geom_point()+
		geom_line()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
		coord_flip()+
		labs(x = "Age (y bp)", y = "Temperature (°C)", title = "Estimated sea surface temperature", color = "Method")+
		facet_nested(~core, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		guides(fill = guide_legend(ncol = 1))

plotA1

sst <- read.csv("data/combined_sst_proxies_separate_columns.csv")
sst %>%
	filter(!is.na(sst_mgca_jonkers_2013) & !is.na(sst_mgca_kozdon_2009)) %>%
	ggplot(., aes(x = sst_mgca_jonkers_2013, y = sst_mgca_kozdon_2009))+
		geom_point()+
		geom_abline(intercept = 0, slope = 1)+
		facet_nested(~core)+
		theme(aspect.ratio = 1)


plotA2 <- cdata %>%
	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202"))) %>%
	ggplot(., aes(x = y_bp, y = mis))+
		geom_point()+
		geom_line()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
		coord_flip()+
		labs(x = "Age (y bp)", y = "Temperature (°C)", title = "Estimated sea surface temperature", color = "Method")+
		facet_nested(~core, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		guides(fill = guide_legend(ncol = 1))


###################################################################
############ Foraminifera data ####################################
###################################################################

foram <- read_tsv("data/combined_foraminifera_geochem.tsv")

plotA <- foram %>%
	select(core, y_bp, polar_planktonic_spp_pct, transitional_planktonic_species_pct, sub_polar_planktonic_species_pct) %>%
	drop_na() %>%
	pivot_longer(cols = -c(core, y_bp),
				 names_to = "observation", values_to = "value") %>%
	mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202")) %>%
	ggplot(., aes(x = y_bp, y = value, fill = observation))+
		geom_area(stat = "identity", position = "stack", color = "black", linewidth = 0.2, alpha = 0.8)+
		coord_flip()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
		facet_nested(~core, scales = "free")+
		scale_fill_manual(values = hcl.colors(n = 5, palette = "Zissou 1")[-2], labels = c("polar_planktonic_spp_pct" = "Polar", "sub_polar_planktonic_species_pct" = "Sub-polar", "transitional_planktonic_species_pct" = "Transitional"))+
		labs(x = "Age (y bp)", y = "Proportion of foraminifera by type (%)", fill = "")+
		guides(color = guide_legend(nrow = 1))

png(file = "plots/proxies/foraminifera_type.png", width = 13, height = 10, units = "in", res = 600)
plotA
dev.off()

egg::ggarrange(plotA1+theme(legend.position = "none"), plotA2, plotA+theme(legend.position = "none"), nrow = 1)

species_names <- c(
  "t_quinqueloba" = "Turborotalita quinqueloba",
  "g_bulloides" = "Globigerina bulloides",
  "g_glutinata" = "Globigerinita glutinata",
  "g_inflata" = "Globorotalia inflata",
  "g_scitula" = "Globorotalia scitula",
  "g_truncatulinoides_dex" = "Globorotalia truncatulinoides",
  "g_uvula" = "Globigerinita uvula",
  "n_incompta" = "Neogloboquadrina incompta",
  "n_pachyderma" = "Neogloboquadrina pachyderma",
  "o_universa" = "Orbulina universa"
)

plotA <- foram %>%
	select(core, y_bp, t_quinqueloba_pct, g_bulloides_pct, g_glutinata_pct, g_glutinata_with_bulla_pct, g_glutinata_without_bulla_pct, g_inflata_pct, g_scitula_pct, g_truncatulinoides_dex, g_uvula_pct, n_incompta_pct, n_pachyderma_pct, o_universa) %>%
	mutate(across(-c(core, y_bp), ~replace_na(., 0))) %>% 
	filter(rowSums(across(-c(core, y_bp))) != 0) %>%
	mutate(g_glutinata_pct = g_glutinata_pct + g_glutinata_with_bulla_pct + g_glutinata_without_bulla_pct) %>%
	select(-g_glutinata_with_bulla_pct, -g_glutinata_without_bulla_pct) %>%
	pivot_longer(cols = -c(core, y_bp),
				 names_to = "observation", values_to = "value") %>%
	mutate(observation = gsub("_pct", "", observation)) %>%
	mutate(observation = recode(observation, !!!species_names)) %>%
	group_by(observation) %>%
	mutate(mean_value = mean(value, na.rm = TRUE)) %>%
	ungroup() %>%
	mutate(observation = fct_reorder(observation, -mean_value)) %>%
	mutate(.core_ybp = paste(core, y_bp, sep = "-")) %>%
	select(-core, -y_bp) %>%
	complete(.core_ybp, observation, fill = list(value = 0)) %>%
	separate(.core_ybp, into = c("core", "y_bp"), sep = "-", extra = "merge", convert = TRUE) %>%
	mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202")) %>%
	ggplot(., aes(x = y_bp, y = value, fill = observation))+
		geom_area(stat = "identity", position = "stack", color = "black", linewidth = 0.2, alpha = 0.8)+
		coord_flip()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
		facet_nested(~core, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		labs(x = "Age (y bp)", y = "Proportion of foraminifera species (%)", fill = "")+
		guides(color = guide_legend(nrow = 1))+
		theme(legend.text = element_text(face = "italic"))
		
png(file = "plots/proxies/foraminifera_sp.png", width = 13, height = 10, units = "in", res = 600)
plotA
dev.off()


################################################
######## DNA library concentration #############
################################################

plotA <- cdata %>%
	filter(!is.na(library_concentration)) %>%
	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
	ggplot(., aes(x = y_bp, y = library_concentration))+
		geom_point()+
		geom_line()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
		coord_flip()+
		labs(x = "Age (y bp)", y = "Library concentration (nM)", color = "Method")+
		facet_nested(~core, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		guides(fill = guide_legend(ncol = 1))

png(file = "plots/proxies/library_dna_concentration.png", width = 13, height = 10, units = "in", res = 600)
plotA
dev.off()

dmg <- read_tsv("/projects/caeg/people/ngm902/apps/repos/rocs/results/microbial/damage/damage-classification-depositional/sample-baselines.tsv")

cdata %>%
	left_join(dmg %>% select(label, sample_A_b_median)) %>%
	filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
	ggplot(., aes(x = library_concentration, y = sample_A_b_median))+
		geom_point()+
		theme(aspect.ratio = 1)+
		scale_x_sqrt()+
		facet_nested(~core)


################################################
######## GeoB25202-specific proxies ############
################################################


geob_proxis <- read_tsv("/projects/caeg/people/ngm902/apps/repos/rocs/data/geob25202_clean_proxies.tsv")

plotA <- geob_proxis %>%
	select(core, y_bp, toc_percent_weight, sum_expandable_clays, sum_illites_micas, sum_montmorillonites_smectites, ratio_greenland_to_iceland) %>%
	mutate(across(everything(), ~replace_na(., 0))) %>%
	pivot_longer(cols = -c(core, y_bp), names_to = "observation", values_to = "value") %>%
	# filter(value != 0) %>%
	ggplot(., aes(x = y_bp, y = value))+
		geom_point()+
		geom_line()+
		coord_flip()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
		facet_nested(~observation, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		labs(x = "Age (y bp)", y = "Value", title = "", fill = "")+
		guides(color = guide_legend(nrow = 1))+
		theme(legend.position = "bottom")

png(file = "plots/proxies/geob25202_additional_proxies.png", width = 13, height = 10, units = "in", res = 600)
plotA
dev.off()


tmpA <- geob_proxis %>%
	select(core, y_bp, toc_percent_weight, sum_expandable_clays, sum_illites_micas, sum_montmorillonites_smectites, ratio_greenland_to_iceland) %>%
	mutate(across(everything(), ~replace_na(., 0))) %>%
	pivot_longer(cols = -c(core, y_bp), names_to = "observation", values_to = "value")

tmpB <- cdata %>%
	filter(core == "GeoB25202_R1") %>%
	mutate(core = gsub("_R1", "", core)) %>%	
	filter(!is.na(temp_method)) %>%
	select(core, y_bp, observation = temp_method, value = temp)

rbind(tmpA, tmpB) %>%
	filter(value != 0) %>%
	ggplot(., aes(x = y_bp, y = value))+
		geom_point()+
		geom_line()+
		coord_flip()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
		facet_nested(~observation, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		labs(x = "Age (y bp)", y = "Value", title = "", fill = "")+
		guides(color = guide_legend(nrow = 1))+
		theme(legend.position = "bottom")


## Correlation between temperature and ratio_greenland_to_iceland??
rbind(tmpA, tmpB) %>%
	filter(observation %in% c("ratio_greenland_to_iceland", "uk37_alkenone")) %>%
	filter(value != 0) %>%
	ggplot(., aes(x = y_bp, y = value))+
		geom_point()+
		geom_line()+
		coord_flip()+
		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
		facet_nested(~observation, scales = "free")+
		scale_fill_manual(values = paired_genus)+
		labs(x = "Age (y bp)", y = "Value", title = "", fill = "")+
		guides(color = guide_legend(nrow = 1))+
		theme(legend.position = "bottom")



library(dplyr)
library(tidyr)
library(fuzzyjoin)

df <- bind_rows(tmpA, tmpB) %>%
  filter(core == "GeoB25202",
         observation %in% c("ratio_greenland_to_iceland", "uk37_alkenone")) %>%
  select(core, y_bp, observation, value)

ratio <- df %>% filter(observation == "ratio_greenland_to_iceland") %>% select(y_bp_ratio = y_bp, ratio = value)
uk37  <- df %>% filter(observation == "uk37_alkenone") %>% select(y_bp_uk37 = y_bp, uk37 = value)

# tolerance in years
tol <- 500

pairs <- difference_left_join(
  uk37, ratio,
  by = c("y_bp_uk37" = "y_bp_ratio"),
  max_dist = tol,
  distance_col = "delta_y"
) %>%
  group_by(y_bp_uk37) %>%
  slice_min(delta_y, n = 1, with_ties = FALSE) %>%  # el más cercano
  ungroup() %>%
  filter(!is.na(ratio))

pairs %>%
  ggplot(aes(x = uk37, y = ratio)) +
	geom_point() +
	geom_smooth(method = "lm", se = TRUE, color = "red") +
	labs(x = "UK37 (°C)", y = "Ratio Greenland/Iceland", title = "Correlation between UK37 and Ratio Greenland/Iceland") +
	ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top")





################################################
############# XRF data #########################
################################################

xrf <- read_tsv("data/combined_xrf_not_normalized.tsv")
# xrf <- read.csv("data/combined_xrf_geochemistry_not_normalized.csv")


# Ratios with smoothing
# windows in years
window_years <- 1000

plotA <- xrf %>%
	filter(y_bp > 0) %>%
	select(-depth_in_core_cm) %>%
	# mutate(across(everything(), ~replace_na(., 0))) %>% 
	select(br, ti, ba, ca, si, p, core, y_bp) %>%
	drop_na() %>%
	mutate(
		ratio_br_ti = br / ti,
		ratio_ba_ti = ba / ti,
		ratio_ca_ti = ca / ti,
		ratio_si_ti = si / ti,
		ratio_p_ti = sqrt(p / ti),
		) %>%
	select(core, y_bp, ratio_br_ti, ratio_ba_ti, ratio_ca_ti, ratio_si_ti, ratio_p_ti) %>%
	pivot_longer(cols = -c(core, y_bp), names_to = "observation", values_to = "value") %>%
	mutate(observation = recode(observation,
	"ratio_br_ti" = "Br/Ti",
	"ratio_ba_ti" = "Ba/Ti",
	"ratio_ca_ti" = "Ca/Ti",
	"ratio_si_ti" = "Si/Ti",
	"ratio_p_ti" = "P/Ti"
	)) %>%
	filter(!is.na(value), is.finite(value)) %>%
	mutate(
		core = fct_relevel(core, "ST5", "ST8", "ST13")
	) %>%
	group_by(core, observation) %>%
	arrange(y_bp, .by_group = TRUE) %>%
	mutate(
		smoothed = {
			# para cada fila, promedia los valores cuya diferencia en y_bp está dentro de la mitad de la ventana
			y_vec <- y_bp
			v_vec <- value
			sapply(y_vec, function(yp) {
				sel <- abs(y_vec - yp) <= (window_years / 2)
				vals <- v_vec[sel]
				vals <- vals[!is.na(vals)]
				if (length(vals) == 0) NA_real_ else mean(vals)
			})
		}
	) %>%
	ungroup() %>%
	ggplot(aes(x = y_bp, color = observation)) +
		geom_line(aes(y = value), alpha = 0.6) +  # línea original con transparencia
		geom_line(aes(y = smoothed), linewidth = 0.2, color = "black") +
		coord_flip() +
		scale_x_reverse(
			labels = scales::label_number(scale_cut = scales::cut_short_scale()),
			breaks = seq(0, 5e06, 2e4)
		) +
		scale_y_continuous(
			labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::breaks_pretty(n = 2)
		) +
		scale_color_manual(values = c("#154360", "#FF5733", "#FFC300", "#1ABC9C", "#8E44AD")) +
		facet_nested(~observation*core, scales = "free") +
		labs(x = "Age (y bp)", y = "Ratio") +
		guides(color = guide_legend(ncol = 1)) +
		theme(legend.position = "none")

png(file = "plots/proxies/XRF_productivity_ratios.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()



# In the late Quaternary Arctic sections, the Al/Si, Fe/Ca, Ti/Al, and Zr/Al ratios were applied to estimate the variability of clastic components [24,25]. 

plotA <- xrf %>%
	select(-depth_in_core_cm) %>%
	# mutate(across(everything(), ~replace_na(., 0))) %>% 
	select(ti, al, fe, ca, zr, rb, sr, core, y_bp) %>%
	drop_na() %>%
	mutate(
		# ratio_al_si  = sqrt(al / si),
		ratio_ti_al  = ti / al,
		ratio_fe_al  = fe / al,
		ratio_fe_ca  = fe / ca,
		ratio_zr_al  = sqrt(zr / al),
		ratio_rb_sr  = rb / sr
	) %>%
	select(core, y_bp, ratio_fe_al, ratio_ti_al, ratio_rb_sr, ratio_fe_ca, ratio_zr_al) %>%
	pivot_longer(cols = -c(core, y_bp), names_to = "observation", values_to = "value") %>%
	mutate(observation = recode(observation,
		"ratio_al_si" = "Al/Si",
		"ratio_ti_al" = "Ti/Al",
		"ratio_fe_al" = "Fe/Al",
		"ratio_rb_sr" = "Rb/Sr",
		"ratio_fe_ca" = "Fe/Ca",
		"ratio_zr_al"  = "Zr/Al"
	)) %>%
	filter(!is.na(value), is.finite(value)) %>%
	mutate(
		observation = factor(observation, levels = c("Al/Si", "Ti/Al", "Fe/Al", "Rb/Sr", "Fe/Ca", "Zr/Al")),
		core = fct_relevel(core, "ST5", "ST8", "ST13")
	) %>%
	group_by(core, observation) %>%
	arrange(y_bp, .by_group = TRUE) %>%
	mutate(
		smoothed = {
			# para cada fila, promedia los valores cuya diferencia en y_bp está dentro de la mitad de la ventana
			y_vec <- y_bp
			v_vec <- value
			sapply(y_vec, function(yp) {
				sel <- abs(y_vec - yp) <= (window_years / 2)
				vals <- v_vec[sel]
				vals <- vals[!is.na(vals)]
				if (length(vals) == 0) NA_real_ else mean(vals)
			})
		}
	) %>%
	ungroup() %>%
	ggplot(aes(x = y_bp, color = observation)) +
		geom_line(aes(y = value), alpha = 0.5) +  # línea original con transparencia
		geom_line(aes(y = smoothed), linewidth = 0.2, color = "black") +
		coord_flip() +
		scale_x_reverse(
			labels = scales::label_number(scale_cut = scales::cut_short_scale()),
			breaks = seq(0, 5e06, 2e4)
		) +
		scale_y_continuous(
			labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::breaks_pretty(n = 3)
		) +
		
		# scale_color_manual(values = wesanderson::wes_palette("Moonrise1")) +
		# scale_color_manual(values = c("#154360", "#FF5733", "#FFC300", "#1ABC9C")) +
		facet_nested(~observation*core, scales = "free") +
		labs(x = "Age (y bp)", y = "Ratio", color = "Ratio") +
		guides(color = guide_legend(ncol = 1)) +
		theme(legend.position = "none")

png(file = "plots/proxies/XRF_terrigenous_ratios.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()



# Manganese, as a redox-sensitive element, was used in the sediments of the Lomonosov Ridge in the Arctic Ocean to determine the conditions of anoxia caused by ice cover [26]. 
# Jakobsson, M.; Løvlie, R.; Al-Hanbali, H.; Arnold, E.; Backman, J.; Mörth, M. Manganese and color cycles in Arctic Ocean sediments constrain Pleistocene chronology. Geology 2000, 28, 23–26. [Google Scholar] [CrossRef]

plotA <- xrf %>%
	select(-depth_in_core_cm) %>%
	# mutate(across(everything(), ~replace_na(., 0))) %>% 
	mutate(
		ratio_mn_fe  = mn / fe,
		mo = sqrt(mo),
		mo_al = sqrt(mo / al),
		mo_ti = mo / ti,
		mn_ti = mn / ti,
		fe_ti = fe / ti,
		si_ti = si / ti
	) %>%
	# select(core, y_bp, mn, mo, mo_al, mn_ti, si_ti, ratio_mn_fe) %>%
	select(core, y_bp, mn_ti, mo_ti, fe_ti, si_ti, ratio_mn_fe) %>%
	pivot_longer(cols = -c(core, y_bp), names_to = "observation", values_to = "value") %>%
	# filter(value != 0) %>%
	mutate(observation = recode(observation,
		"mn_ti" = "Mn/Ti",
		"si_ti" = "Si/Ti",
		"mo_ti" = "Mo/Ti",
		"fe_ti" = "Fe/Ti",
		# "mo_al" = "Mo/Al",
		"ratio_mn_fe" = "Mn/Fe",
		# "mn" = "Mn",
		# "mo" = "Mo",
	)) %>%
	filter(!is.na(value), is.finite(value)) %>%
	mutate(core = fct_relevel(core, "ST5", "ST8", "ST13")) %>%
	group_by(core, observation) %>%
	arrange(y_bp, .by_group = TRUE) %>%
	mutate(
		smoothed = {
			y_vec <- y_bp
			v_vec <- value
			sapply(y_vec, function(yp) {
				sel <- abs(y_vec - yp) <= (window_years / 2)
				vals <- v_vec[sel]
				vals <- vals[!is.na(vals)]
				if (length(vals) == 0) NA_real_ else mean(vals)
			})
		}
	) %>%
	ungroup() %>%
	ggplot(aes(x = y_bp, color = observation)) +
		geom_line(aes(y = value), alpha = 0.5) + 
		geom_line(aes(y = smoothed), linewidth = 0.2, color = "black") +
		coord_flip() +
		scale_x_reverse(
			labels = scales::label_number(scale_cut = scales::cut_short_scale()),
			breaks = seq(0, 5e06, 2e4)
		) +
		scale_y_continuous(
			labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::breaks_pretty(n = 3)
		) +
		facet_nested(~observation*core, scales = "free") +
		labs(x = "Age (y bp)", y = "Ratio", color = "Ratio") +
		guides(color = guide_legend(ncol = 1)) +
		theme(legend.position = "none")

png(file = "plots/proxies/XRF_redox_ratios.png", width = 15, height = 10, units = "in", res = 600)
plotA
dev.off()
















########################################################################
######### Function to perform PCA on XRF data and plot results #########
########################################################################

library(dplyr)
library(compositions)
library(ggplot2)
library(factoextra)
library(ggrepel)
library(egg)
library(tibble)
library(scales)
library(grid)

vars_xrf <- c("al","fe","k","rb","ti","zr","ca","sr","br","ba","mn","s","mo")

replace_small <- function(x) {
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[is.na(x)] <- min_pos / 2
  x[x == 0]   <- min_pos / 2
  x
}


make_pca_plot <- function(core_name, xrf, outdir = "plots/proxies", arrow_scale = 15) {
  
#   dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  xrf_clr <- xrf %>%
    filter(core == core_name) %>%
    select(y_bp, all_of(vars_xrf)) %>%
	drop_na()
    # mutate(across(all_of(vars_xrf), replace_small))
  
  # Comprobar que hay datos
  if (nrow(xrf_clr) < 3) {
    message("Saltando ", core_name, ": no hay suficientes datos.")
    return(NULL)
  }
  
  mat <- xrf_clr %>%
    select(all_of(vars_xrf)) %>%
    as.matrix()
  
  mat_clr <- clr(mat)
  
  pca_clr <- prcomp(mat_clr, center = TRUE, scale. = TRUE)
  
  scores_clr <- as.data.frame(pca_clr$x) %>%
    bind_cols(xrf_clr %>% select(y_bp))
  
  loadings_clr <- as.data.frame(pca_clr$rotation) %>%
    rownames_to_column("element")
  
  loadings_plot <- loadings_clr %>%
    mutate(
      PC1 = PC1 * arrow_scale,
      PC2 = PC2 * arrow_scale
    )
  
  plotA <- ggplot() +
    geom_point(
      data = scores_clr,
      aes(PC1, PC2, color = y_bp),
      size = 1.5,
      alpha = 0.7
    ) +
    geom_segment(
      data = loadings_plot,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text_repel(
      data = loadings_plot,
      aes(x = PC1, y = PC2, label = element),
      size = 3
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = paste("PCA XRF -", core_name),
      x = "PC1",
      y = "PC2",
      color = "Age (yr BP)"
    )
  
  plotB <- scores_clr %>%
	ggplot(aes(y = PC1, x = y_bp, color = y_bp)) +
	geom_line() +
	coord_flip() +
	scale_color_viridis_c(direction = -1) +
	labs(y = "PC1", x = "Age (y bp)", color = "Age") +
	theme_bw() +
	theme(legend.position = "none") +
	scale_x_reverse(
	breaks = function(x) pretty(x, n = 10),  # ← clave
	labels = scales::label_number(scale_cut = scales::cut_short_scale())
	)
  
 plotC <- scores_clr %>%
	ggplot(aes(y = PC2, x = y_bp, color = y_bp)) +
	geom_line() +
	coord_flip() +
	scale_color_viridis_c(
		direction = -1,
		labels = scales::label_number(scale_cut = scales::cut_short_scale())
	) +
	labs(y = "PC2", x = NULL, color = "Age") +  # ← sin título
	theme_bw() +
	theme(
	# legend.position = "none",
	axis.text.y = element_blank(),   # ← sin números
	axis.ticks.y = element_blank()   # ← sin ticks
	) +
	scale_x_reverse()
  
  combined_plot <- egg::ggarrange(
    plotA, plotB, plotC,
    nrow = 1,
    widths = c(2, 1, 1)
  )
  
  outfile <- file.path(outdir, paste0("PCA_XRF_", core_name, ".png"))
  


  png(filename = outfile, width = 13, height = 6, units = "in", res = 600)
  print(combined_plot)
  dev.off()
  
  message("Guardado: ", outfile)
  
  invisible(list(
    core = core_name,
    pca = pca_clr,
    scores = scores_clr,
    loadings = loadings_clr
  ))
}

# Cores a exportar
cores_to_plot <- c("ST5", "ST8", "ST13")

# Ejecutar para cada core
results <- lapply(cores_to_plot, make_pca_plot, xrf = xrf)
names(results) <- cores_to_plot












########################################################################
######### Function to perform PCA on XRF data and plot results #########
########################################################################

library(dplyr)
library(compositions)
library(ggplot2)
library(factoextra)
library(ggrepel)
library(egg)
library(tibble)
library(scales)
library(grid)


	
vars_xrf <- c("ratio_ti_al", "ratio_fe_al", "ratio_rb_sr", "ratio_fe_ca", "ratio_zr_al")

replace_small <- function(x) {
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[is.na(x)] <- min_pos / 2
  x[x == 0]   <- min_pos / 2
  x
}

make_pca_plot <- function(core_name, xrf, outdir = "plots/proxies", arrow_scale = 15) {
  
  xrf_clr <- xrf %>%
    filter(core == core_name) %>%
	# mutate(across(everything(), ~replace_na(., 0))) %>% 
	select(y_bp, core, ti, al, fe, ca, zr, rb, sr) %>%
	drop_na() %>%
	mutate(
		ratio_ti_al  = ti / al,
		ratio_fe_al  = fe / al,
		ratio_rb_sr  = ifelse(core != "GeoB25202", (rb / sr), rb_sr),
		ratio_fe_ca  = fe / ca,
		ratio_zr_al  = sqrt(zr / al)
	) %>%
    select(y_bp, all_of(vars_xrf)) %>%
    mutate(across(all_of(vars_xrf), replace_small))
  
  # Comprobar que hay datos
  if (nrow(xrf_clr) < 3) {
    message("Saltando ", core_name, ": no hay suficientes datos.")
    return(NULL)
  }
  
  mat <- xrf_clr %>%
    select(all_of(vars_xrf)) %>%
    as.matrix()
  
  mat_clr <- clr(mat)
  
  pca_clr <- prcomp(mat_clr, center = TRUE, scale. = TRUE)
  
  scores_clr <- as.data.frame(pca_clr$x) %>%
    bind_cols(xrf_clr %>% select(y_bp))
  
  loadings_clr <- as.data.frame(pca_clr$rotation) %>%
    rownames_to_column("element")
  
  loadings_plot <- loadings_clr %>%
    mutate(
      PC1 = PC1 * arrow_scale,
      PC2 = PC2 * arrow_scale
    ) %>%
	mutate(element = recode(element,
		"ratio_ti_al" = "Ti/Al",
		"ratio_fe_al" = "Fe/Al",
		"ratio_rb_sr" = "Rb/Sr",
		"ratio_fe_ca" = "Fe/Ca",
		"ratio_zr_al"  = "Zr/Al"
	))
  
  plotA <- ggplot() +
    geom_point(
      data = scores_clr,
      aes(PC1, PC2, color = y_bp),
      size = 1.5,
      alpha = 0.7
    ) +
    geom_segment(
      data = loadings_plot,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text_repel(
      data = loadings_plot,
      aes(x = PC1, y = PC2, label = element),
      size = 3
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = paste("PCA XRF -", core_name),
      x = "PC1",
      y = "PC2",
      color = "Age (yr BP)"
    )
  
  plotB <- scores_clr %>%
	ggplot(aes(y = PC1, x = y_bp, color = y_bp)) +
	geom_line() +
	coord_flip() +
	scale_color_viridis_c(direction = -1) +
	labs(y = "PC1", x = "Age (y bp)", color = "Age") +
	theme_bw() +
	theme(legend.position = "none") +
	scale_x_reverse(
	breaks = function(x) pretty(x, n = 10),  # ← clave
	labels = scales::label_number(scale_cut = scales::cut_short_scale())
	)

 plotC <- scores_clr %>%
	ggplot(aes(y = PC2, x = y_bp, color = y_bp)) +
	geom_line() +
	coord_flip() +
	scale_color_viridis_c(
		direction = -1,
		labels = scales::label_number(scale_cut = scales::cut_short_scale())
	) +
	labs(y = "PC2", x = NULL, color = "Age") +  # ← sin título
	theme_bw() +
	theme(
	# legend.position = "none",
	axis.text.y = element_blank(),   # ← sin números
	axis.ticks.y = element_blank()   # ← sin ticks
	) +
	scale_x_reverse()
  
  combined_plot <- egg::ggarrange(
    plotA, plotB, plotC,
    nrow = 1,
    widths = c(2, 1, 1)
  )
  
  outfile <- file.path(outdir, paste0("PCA_XRF_", core_name, "-TERRIGENOUS.png"))
  


  png(filename = outfile, width = 13, height = 6, units = "in", res = 600)
  print(combined_plot)
  dev.off()
  
  message("Guardado: ", outfile)
  
  invisible(list(
    core = core_name,
    pca = pca_clr,
    scores = scores_clr,
    loadings = loadings_clr
  ))
}

# Cores a exportar
cores_to_plot <- c("ST5", "ST8", "ST13")

# Ejecutar para cada core
results <- lapply(cores_to_plot, make_pca_plot, xrf = xrf)
names(results) <- cores_to_plot










########################################################################
######### Function to perform PCA on XRF data and plot results #########
########################################################################

library(dplyr)
library(compositions)
library(ggplot2)
library(factoextra)
library(ggrepel)
library(egg)
library(tibble)
library(scales)
library(grid)
		
vars_xrf <- c("ratio_ba_ti", "ratio_br_ti", "ratio_ca_ti", "ratio_p_ti", "ratio_si_ti")

replace_small <- function(x) {
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[is.na(x)] <- min_pos / 2
  x[x == 0]   <- min_pos / 2
  x
}

# # aplicar smoothing por medias móviles en una ventana temporal (años)
window_years_local <- if (exists("window_years")) window_years else 1000

# xrf %>%
# 	filter(core == "ST13") %>%
# 	select(y_bp, p, ti) %>%
# 	arrange(y_bp) %>%
# 	mutate(
# 		# media móvil basada en ventana temporal alrededor de cada punto
# 		smoothed_p = {
# 			y_vec <- y_bp
# 			v_vec <- p
# 			sapply(y_vec, function(yp) {
# 				sel <- abs(y_vec - yp) <= (window_years_local / 2)
# 				vals <- v_vec[sel]
# 				vals <- vals[!is.na(vals)]
# 				if (length(vals) == 0) NA_real_ else mean(vals)
# 			})
# 		},
# 		# evitar ceros si es necesario
# 		smoothed_p = ifelse(is.na(smoothed_p) & !is.na(p), p, smoothed_p),
# 		ratio = smoothed_p / ti
# 	) %>%
# 	mutate(p = ifelse(p == 0, 1e-6, p)) %>%
# 	ggplot(aes(x = y_bp)) +
# 		geom_point(aes(y = p * 100), color = "lightblue", alpha = 0.5) +   # p original (transparente)
# 		geom_point(aes(y = smoothed_p * 100), color = "blue") +            # p suavizada
# 		geom_point(aes(y = ti), color = "red") +
# 		geom_line(aes(y = ratio * 500), color = "black") +
# 		geom_point(aes(y = ratio * 500), color = "black") +
# 		scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))


# xrf_clr %>%
# 	# pull(ratio_p_ti) %>% summary()
# 	ggplot(aes(x = y_bp)) +
# 		geom_point(aes(y = p * 100), color = "lightblue", alpha = 0.5) +   # p original (transparente)
# 		geom_point(aes(y = smoothed_p * 100), color = "blue") +            # p suavizada
# 		geom_point(aes(y = ti), color = "red") +
# 		geom_line(aes(y = ratio_p_ti * 500), color = "black") +
# 		geom_point(aes(y = ratio_p_ti * 500), color = "black") +
# 		scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))



core_name <- "ST13"

make_pca_plot <- function(core_name, xrf, outdir = "plots/proxies", arrow_scale = 10) {
  
  window_years_local <- if (exists("window_years")) window_years else 10000

  xrf_clr <- xrf %>%
	filter(core == core_name) %>%
	select(y_bp, core, ti, ba, br, ca, si, p) %>%
	drop_na() %>%
	arrange(y_bp) %>%
	# calcular p suavizada usando la misma ventana temporal
	mutate(
	  smoothed_p = {
		y_vec <- y_bp
		v_vec <- p
		sapply(y_vec, function(yp) {
		  sel <- abs(y_vec - yp) <= (window_years_local / 2)
		  vals <- v_vec[sel]
		  vals <- vals[!is.na(vals)]
		  if (length(vals) == 0) NA_real_ else mean(vals)
		})
	  },
	  # si no hay smoothing disponible, conservar p original cuando exista
	  smoothed_p = ifelse(is.na(smoothed_p) & !is.na(p), p, smoothed_p)
	) %>%
	mutate(
	  ratio_ba_ti = ba / ti,
	  ratio_br_ti = br / ti,
	  ratio_ca_ti = ca / ti,
	  ratio_si_ti = si / ti,
	  # usar p suavizada para el ratio P/Ti
	  ratio_p_ti  = smoothed_p / ti
	) %>%
	select(y_bp, all_of(vars_xrf))  %>%
	mutate(across(all_of(vars_xrf), replace_small))
  
  # Comprobar que hay datos
  if (nrow(xrf_clr) < 3) {
    message("Saltando ", core_name, ": no hay suficientes datos.")
    return(NULL)
  }
  
  mat <- xrf_clr %>%
    select(all_of(vars_xrf)) %>%
    as.matrix()
  
  mat_clr <- clr(mat)
  
  pca_clr <- prcomp(mat_clr, center = TRUE, scale. = TRUE)
  
  scores_clr <- as.data.frame(pca_clr$x) %>%
    bind_cols(xrf_clr %>% select(y_bp))
  
  loadings_clr <- as.data.frame(pca_clr$rotation) %>%
    rownames_to_column("element")
  
  loadings_plot <- loadings_clr %>%
    mutate(
      PC1 = PC1 * arrow_scale,
      PC2 = PC2 * arrow_scale
    ) %>%
	mutate(element = recode(element,
		"ratio_br_ti" = "Br/Ti",
		"ratio_ba_ti" = "Ba/Ti",
		"ratio_ca_ti" = "Ca/Ti",
		"ratio_si_ti" = "Si/Ti",
		"ratio_p_ti" = "P/Ti"
	))
  
scores_clr %>%
  count(y_bp) %>%
  filter(n > 1)

  plotA <- ggplot() +
    geom_point(
      data = scores_clr,
      aes(PC1, PC2, color = y_bp),
      size = 1.5,
      alpha = 0.7
    ) +
    geom_segment(
      data = loadings_plot,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text_repel(
      data = loadings_plot,
      aes(x = PC1, y = PC2, label = element),
      size = 3
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = paste("PCA XRF -", core_name),
      x = "PC1",
      y = "PC2",
      color = "Age (yr BP)"
    )
  
  plotB <- scores_clr %>%
	ggplot(aes(y = PC1, x = y_bp, color = y_bp)) +
	geom_line() +
	coord_flip() +
	scale_color_viridis_c(direction = -1) +
	labs(y = "PC1", x = "Age (y bp)", color = "Age") +
	theme_bw() +
	theme(legend.position = "none") +
	scale_x_reverse(
	breaks = function(x) pretty(x, n = 10),  # ← clave
	labels = scales::label_number(scale_cut = scales::cut_short_scale())
	)

 plotC <- scores_clr %>%
	ggplot(aes(y = PC2, x = y_bp, color = y_bp)) +
	geom_line() +
	coord_flip() +
	scale_color_viridis_c(
		direction = -1,
		labels = scales::label_number(scale_cut = scales::cut_short_scale())
	) +
	labs(y = "PC2", x = NULL, color = "Age") +  # ← sin título
	theme_bw() +
	theme(
	# legend.position = "none",
	axis.text.y = element_blank(),   # ← sin números
	axis.ticks.y = element_blank()   # ← sin ticks
	) +
	scale_x_reverse()
  
  combined_plot <- egg::ggarrange(
    plotA, plotB, plotC,
    nrow = 1,
    widths = c(2, 1, 1)
  )
  
  outfile <- file.path(outdir, paste0("PCA_XRF_", core_name, "-PRODDUCTIVITY.png"))
  
  png(filename = outfile, width = 13, height = 6, units = "in", res = 600)
  print(combined_plot)
  dev.off()
  
  message("Guardado: ", outfile)
  
  invisible(list(
    core = core_name,
    pca = pca_clr,
    scores = scores_clr,
    loadings = loadings_clr
  ))
}

# Cores a exportar
cores_to_plot <- c("ST5", "ST8", "ST13", "GeoB25202")

# Ejecutar para cada core
results <- lapply(cores_to_plot, make_pca_plot, xrf = xrf)
names(results) <- cores_to_plot








