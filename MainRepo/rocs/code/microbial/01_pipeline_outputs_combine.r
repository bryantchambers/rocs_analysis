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

httpgd::hgd()
httpgd::hgd_browse()

# Load labels data
cdata <- read_tsv(file = "data/metadata_v5.tsv")


# Load lca
samples <- cdata$label |> unique()

tax_data <- read_tsv("taxonomic-profiling/microbial/taxonomic-profiling-reassign-summary/standard/k1000/94/tp-lca.summary.tsv.gz") %>%
    filter(rank %in% c("subspecies", "species", "genus")) %>% 
    separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "specie", "subspecie"), sep = ";") %>%
    filter(domain != "d__Eukaryota") %>%
    inner_join(cdata)


# Let's get the damage at each node
dmg <- read_tsv("taxonomic-profiling/microbial/taxonomic-profiling-reassign-damage-lca/standard/k1000/94/tp-damage.tsv.gz") %>%
    filter(label %in% samples) |>
    inner_join(tax_data |> select(label, taxid, name, rank))


# For subspecies
# Let's identify all those taxa with a bad fit
# We set a damage
dmg_ssp <- dmg |>
    filter(rank == "subspecies") |>
    filter(nreads >= 100) |>
    rename(tax_name = taxid, n_reads = nreads)

dat_ssp <- dmg_fwd_CCC(dmg_ssp, samples, ci = "asymptotic", nperm = 100, nproc = 24)

dat1ssp <- dat_ssp |>
    ungroup() |>
    inner_join(dmg_ssp) |>
    mutate(fit = ifelse(rho_c >= 0.75 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |> ## I LOWERED THE rho_c FROM 0.85 TO 0.75
    mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))

dat1ssp %>%
    write_tsv(., "./results/microbial/damage/dmg-subspecies.tsv.gz")

# dat1ssp <- read_tsv("./results/microbial/damage/dmg-subspecies.tsv.gz")

dat1ssp %>%
    select(label, tax_name, n_reads, A_b, q_b, c_b, phi_b, Zfit, fit, rho_c, C_b, rho_c_perm_pval, q_CI_h, c_CI_l) |>
    write_tsv("./results/microbial/damage/dmg-subspecies-lite.tsv.gz")


#### Additionaly to the fit to the expected damage pattern, should we set a damage threshold? ####

## At the suspecies level, the lca and mapping abundance should be equivalent
# LCA to get the taxonomic path and node-level abundance
lca_ssp <- read_tsv( "taxonomic-profiling/microbial/taxonomic-profiling-reassign-summary/standard/k1000/94/tp-lca.summary.tsv.gz") %>%
    filter(rank == "subspecies") |>
    separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"), sep = ";") %>%
    filter(domain != "d__Eukaryota") %>%
    select(-n_reads, -rank)

# Mapping stats from filterBAM
mapping_stats <- read_tsv("taxonomic-profiling/microbial/taxonomic-profiling-reassign-summary/standard/k1000/94/tp-mapping-filtered.summary.tsv.gz") |>
    mutate(name = paste0("S__", reference)) |>
    rename(n_reads_ref = n_reads) 

# Combine LCA and mapping stats
tax_data_ssp <- lca_ssp |>
    inner_join(mapping_stats)

# Add damage and LCA reads. If the fit is bad, the damage will be 0
dmg_ssp <- read_tsv("./results/microbial/damage/dmg-subspecies-lite.tsv.gz") |>
    rename(taxid = tax_name) #|>
    # mutate(damage = ifelse(fit == "bad", 0, A_b))

# Combine damage and taxonomic assignments
tax_data_dmg <- tax_data_ssp %>% 
    inner_join(dmg_ssp) %>%
    left_join(cdata)


tax_data_dmg %>%
    write_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")



# tax_data_dmg




# tax_data_dmg %>%
#     filter()


# tax_data_dmg %>%
#     filter(y_bp > 100000) %>%
#     filter(core == "GeoB25202_R1") %>%
#     ggplot(., aes(x = read_length_mean, y = A_b, color = Zfit >= 2, size = n_reads))+
#         geom_point()+
#         geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
#         facet_wrap(~core)


# df <- tax_data_dmg %>%
#   filter(y_bp > 100000, core == "GeoB25202_R1")

# plotly::plot_ly(
#   df,
#   x = ~read_length_mean,
#   y = ~A_b,
#   z = ~rho_c,
# #   type = "scatter3d",
# #   mode = "markers",
#   color = ~(Zfit >= 2)
# #   marker = list(size = ~pmax(2, sqrt(n_reads)))
# ) %>%
#   plotly::layout(
#     scene = list(
#       xaxis = list(title = "read_length_mean"),
#       yaxis = list(title = "A_b"),
#       zaxis = list(title = "rho_c")
#     )
#   )



# plotA <- tax_data_dmg %>%
#     filter(y_bp > 10000) %>%
#     mutate(tax = paste(domain, phylum, class, order, family, genus, species, subspecies)) %>%
#     group_by(tax) %>%
#     mutate(median_dmg = median(A_b)) %>%
#     ungroup() %>%
#     mutate(tax = fct_reorder(tax, median_dmg)) %>%
#     mutate(tax_position = as.integer(tax)) %>% 
#     mutate(zfit_significance = ifelse(Zfit >= 2, "Significant", "Non-significant")) %>%
#     ggplot(., aes(x = tax_position, y = A_b, fill = zfit_significance))+
#         geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
#         # geom_point(alpha = 0.2, aes(color = zfit_significance))+
#         geom_point(alpha = 0.2, shape = 21, color = "black")+
#         labs(x = "Subspecies sorted by median damage", y = "Damage", fill = "Damage significance (Zfit > 2)")+
#         theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
#         # geom_boxplot(aes(group = tax_position))+
#         # scale_y_sqrt(limits = c(0, 0.75), breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8))

# plotA



# ## Damaged and non-damaged determination

# tax_data_euk <- read_tsv("taxonomic-profiling/microbial/taxonomic-profiling-reassign-summary/standard/k1000/94/tp-lca.summary.tsv.gz") %>%
#     filter(rank %in% c("subspecies", "species", "genus")) %>% 
#     separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "specie", "subspecie"), sep = ";") %>%
#     filter(domain == "d__Eukaryota") %>%
#     inner_join(cdata)

# tax_data_euk %>% head()
# # Let's get the damage at each node
# dmg_euk <- read_tsv("taxonomic-profiling/microbial/taxonomic-profiling-reassign-damage-lca/standard/k1000/94/tp-damage.tsv.gz") %>%
#     filter(label %in% samples) |>
#     inner_join(tax_data_euk |> select(label, taxid, name, rank, domain, lineage, kingdom, phylum, class, order, family, genus, specie, subspecie))


# dmg_ssp_euk <- dmg_euk |>
#     filter(rank == "subspecies") |>
#     filter(nreads >= 100) |>
#     rename(tax_name = taxid, n_reads = nreads)

# dmg_ssp_euk %>% head() %>% View()
# dmg_ssp_euk %>% 
#     filter(!grepl(".mito|.plas", name)) %>% 
#     filter(Zfit >= 2) %>%
#     left_join(cdata) %>%
#     mutate(Zfit_significance = ifelse(Zfit >= 2, "Significant", "Non-significant")) %>%
#     ggplot(., aes(x = y_bp, y = A_b, color = Zfit_significance))+
#         geom_point()+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         coord_flip()+
#         facet_nested(~core)



# # Get the LV values
# # From https://mgimond.github.io/ES218/Week08b.html
# lsum <- function(x, l = 6, all = TRUE) {
#   # Limit max letter summaries to 9
#   if (l > 10) {
#     print("Limit level summary to 9")
#     return()
#   }
#   # letter summary labels
#   let <- c("M", "H", "E", "D", "C", "B", "A", "Z", "Y", "X")
#   # Remove missing values
#   x <- na.omit(x)
#   # Sort values
#   x <- sort(x)
#   # Find depths from each end
#   n <- length(x)
#   Lrnk <- vector()
#   Mrnk <- vector()
#   Rrnk <- vector()
#   Lrnk[1] <- n
#   Mrnk[1] <- n
#   Rrnk[1] <- n
#   i <- 1
#   while ((i <= l) & (Lrnk[i] > 1)) {
#     i <- i + 1
#     Lrnk[i] <- floor(Lrnk[i - 1] + 1) / 2
#     Mrnk[i] <- floor(Lrnk[i])
#     Rrnk[i] <- floor(Lrnk[i] + 0.5)
#   }
#   # Get final set of letters
#   val <- factor(let[1:length(Lrnk[-1])], levels = let[1:length(Lrnk[-1])])
#   # Find the summary values
#   LO <- (x[Mrnk[-1]] + x[Rrnk[-1]]) / 2
#   HI <- (x[n - Mrnk[-1] + 1] + x[n - Rrnk[-1] + 1]) / 2
#   MD <- (LO + HI) / 2
#   SP <- HI - LO
#   # Generate output
#   if (all == TRUE) {
#     out <- data.frame(
#       letter = val, depth = Lrnk[-1], lower = LO,
#       mid = MD, upper = HI, spread = SP
#     )
#   } else {
#     out <- data.frame(letter = val, mid = MD)
#   }
#   return(out)
# }



# ab_limits_df <- dmg_ssp_euk %>% 
#   filter(!grepl(".mito|.plas", name)) %>%
#   filter(Zfit >= 2) %>%
#   left_join(cdata, by = join_by(label)) %>%
#   group_by(label) %>%
#   summarise(
#     A_b_limit = {
#       vals <- na.omit(A_b)

#       # lsum() likely needs more than 1 observation
#       if (length(vals) < 2) {
#         NA_real_
#       } else {
#         tmp <- lsum(vals, l = 8)
#         out <- tmp$lower[tmp$letter == "C"]
#         if (!length(out)) NA_real_ else as.numeric(out[[1]])
#       }
#     },
#     .groups = "drop"
#   )


# ab_limits_df %>%
#     left_join(cdata) %>%
#     ggplot(., aes(x = y_bp, y = A_b_limit))+
#         geom_point()+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         coord_flip()+
#         facet_nested(~core)


# ab_limits_df <- dmg_ssp_euk %>% 
#     filter(!grepl(".mito|.plas", name)) %>%
#     filter(Zfit >= 2) %>%
#     left_join(cdata) %>%
#     group_by(label) %>%
#     summarise(
#         A_b_limit = {
#             vals <- na.omit(A_b)
#             if(length(vals) == 0) {
#                 NA_real_
#             } else {
#                 tmp <- lsum(vals, l = 8)
#                 if(!"C" %in% tmp$letter) {
#                     NA_real_
#                 } else {
#                     tmp %>% filter(letter == "C") %>% pull(lower)
#                 }
#             }
#         },
#         .groups = "drop"
#     )

# head(ab_limits_df)



#     head()
    

# dmg_ssp_euk %>% 
#     filter(!grepl(".mito|.plas", name)) %>%
#     group_by(domain, lineage, kingdom, phylum, class, order, family, genus) %>%
#     summarise(n_reads = sum(n_reads), .groups = "drop")  %>%
#     arrange(desc(n_reads)) %>%
#     # select(-n_reads) %>%
#     print(n = Inf)



# plotA <- tax_data_dmg %>%
#     mutate(tax = paste(domain, phylum, class, order, family, genus)) %>%
#     group_by(tax) %>%
#     mutate(median_dmg = median(A_b)) %>%
#     ungroup() %>%
#     mutate(tax = fct_reorder(tax, median_dmg)) %>%
#     mutate(tax_position = as.integer(tax)) %>% 
#     mutate(zfit_significance = ifelse(Zfit >= 2, "Significant", "Non-significant")) %>%
#     ggplot(., aes(x = tax_position, y = A_b, fill = zfit_significance))+
#         geom_hline(yintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.5)+
#         # geom_point(alpha = 0.2, aes(color = zfit_significance))+
#         geom_point(alpha = 0.2, shape = 21, color = "black")+
#         labs(x = "Subspecies sorted by median damage", y = "Damage", fill = "Damage significance (Zfit > 2)")+
#         theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
#         # geom_boxplot(aes(group = tax_position))+
#         # scale_y_sqrt(limits = c(0, 0.75), breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8))

# plotA












# # tax_data_ssp <- lca_ssp |>
# #     select(-nreads, -abundance) |>
# #     inner_join(mapping_stats)

# # Add dmage. If the fit is bad, the damage will be 0
# dmg_ssp <- read_tsv("./results/microbial/damage/dmg-subspecies-lite.tsv.gz") |>
#     rename(taxid = tax_name) |>
#     mutate(damage = ifelse(fit == "bad", 0, A_b))

# # Combine damage and taxonomic assignments
# tax_data_dmg <- tax_data_ssp %>% 
#     inner_join(dmg_ssp)

# tax_data_dmg %>% names()

# # IF WE SET A DAMAGE THRESHOLD WE SHOULD INCLUDE IT HERE

# # We defined a taxa is damaged if the damage is above the 0.1
# tax_data_dmg <- tax_data_dmg |>
#     ungroup() |>
#     left_join(cdata) %>%
#     mutate(is_dmg = ifelse(damage >= 0.1 | (damage >= 0.05 & y_bp < 10000), "Damaged", "Non-damaged")) |>
#     # mutate(is_dmg = ifelse(damage >= 0.05 & y_bp < 10000, "Damaged", is_dmg)) |>
#     mutate(is_dmg = ifelse(Zfit < 2, "Non-damaged", is_dmg)) %>%
#     mutate(is_dmg = ifelse(fit == "bad", "Non-damaged", is_dmg))

# tax_data_dmg |>
#     write_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")





# tax_data_dmg <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")


# tax_data_dmg %>%
#     # filter(breadth >= 0.01) %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     filter(y_bp > 350000) %>%
#     filter(label %in% unique(label)[1:10]) %>%
#     ggplot(., aes(x = read_ani_mean, y = A_b, color = is_dmg))+
#         geom_point()+
#         facet_wrap(~label)
    




# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = y_bp, y = A_b, color = is_dmg))+
#         geom_point(alpha = 0.01)+
#         geom_line(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), .groups = "drop"), 
#             aes(y = mean_A_b), linewidth = 1)+
#         geom_point(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), .groups = "drop"), 
#             aes(y = mean_A_b, fill = is_dmg), shape = 21, color = "black", size = 2)+
#         labs(y = "Damage (A)", x = "Age (y bp)", fill = "References classified as", color = "References classified as")+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         scale_y_continuous(breaks = c(0, .2, .4, .6, .8))+
#         facet_nested(~core) +
#         coord_flip()+
#         theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 6/1)

# png(file = "/projects/caeg/people/ngm902/apps/repos/rocs_marine/plots/microbial/A_b_and_damage_s.png", width = 8, height = 9, units = "in", res = 600)
# plotA
# dev.off()


# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = A_b, y = y_bp, color = is_dmg))+
#         geom_point(size = 1, alpha = 0.01)+
#         geom_point(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), sd_A_b = sd(A_b), .groups = "drop"), 
#             aes(x = mean_A_b, fill = is_dmg), size = 3, shape = 21, color = "black") +
#         geom_errorbarh(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), n = sum(!is.na(A_b)), sd_A_b = sd(A_b)/sqrt(n), .groups = "drop"), 
#             aes(x = mean_A_b, xmin = mean_A_b - sd_A_b, xmax = mean_A_b + sd_A_b), color = "black", height = 5, alpha = 1) +
#         labs(x = "Damage (A)", y = "Age (y bp)", fill = "References classified as", color = "References classified as")+
#         scale_y_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         facet_nested(~core) +
#         theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 5/1)

# png(file = "/projects/caeg/people/ngm902/apps/repos/rocs_marine/plots/microbial/A_b_and_damage_s.png", width = 8, height = 9, units = "in", res = 600)
# plotA
# dev.off()

# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = A_b, y = y_bp, color = is_dmg))+
#         geom_point(size = 1, alpha = 0.6)+
#         geom_point(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), sd_A_b = sd(A_b), .groups = "drop"), 
#             aes(x = mean_A_b, fill = is_dmg), size = 3, shape = 21, color = "black") +
#         geom_errorbarh(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), n = sum(!is.na(A_b)), sd_A_b = sd(A_b)/sqrt(n), .groups = "drop"), 
#             aes(x = mean_A_b, xmin = mean_A_b - sd_A_b, xmax = mean_A_b + sd_A_b), color = "black", height = 5, alpha = 1) +
#         labs(x = "Damage (A_b)", y = "Age (y bp)", fill = "References classified as", color = "References classified as")+
#         scale_y_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         facet_nested(~core) +
#         theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 5/1)

# png(file = "/projects/caeg/people/ngm902/apps/repos/rocs_marine/plots/microbial/A_b_and_damage.png", width = 10, height = 10, units = "in", res = 600)
# plotA
# dev.off()


# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = A_b, y = y_bp, color = is_dmg))+
#         geom_point(size = 1, alpha = 0.1)+
#         # geom_point(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), sd_A_b = sd(A_b), .groups = "drop"), 
#         #     aes(x = mean_A_b, fill = is_dmg), size = 3, shape = 21, color = "black") +
#         # geom_errorbarh(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), n = sum(!is.na(A_b)), sd_A_b = sd(A_b)/sqrt(n), .groups = "drop"), 
#         #     aes(x = mean_A_b, xmin = mean_A_b - sd_A_b, xmax = mean_A_b + sd_A_b), color = "black", height = 5, alpha = 1) +
#         geom_line(data = . %>% group_by(core, y_bp, is_dmg) %>% summarise(mean_A_b = median(A_b), sd_A_b = sd(A_b), .groups  = "drop"), 
#             aes(x = mean_A_b, color = is_dmg)) +
#         labs(x = "Damage (A_b)", y = "Age (y bp)", fill = "References classified as", color = "References classified as")+
#         scale_y_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         facet_nested(~core) +
#         theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 5/1)

# png(file = "/projects/caeg/people/ngm902/apps/repos/rocs_marine/plots/microbial/A_b_and_damage.png", width = 10, height = 10, units = "in", res = 600)
# plotA
# dev.off()



# # Weighted number of species damaged vs non-damaged
# nsp_fit <- tax_data_dmg %>%
#   group_by(label) %>%
#   count() %>%
#   ungroup()

# plotA <- tax_data_dmg %>%
#   group_by(label, is_dmg) |>
#   summarise(n_reads = sum(n_reads), .groups = "drop") %>%
#   group_by(label) %>%
#   mutate(relative_n_reads = n_reads / sum(n_reads)) %>%
#   ungroup() %>% 
#   left_join(nsp_fit) %>%
#   mutate(proportion = n*relative_n_reads) %>%
#   select(label, is_dmg, proportion) %>%
#   complete(label, is_dmg, fill = list(proportion = 0)) %>%
#   inner_join(cdata) |> 
#   mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#   mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#   ggplot(aes(x = y_bp, y = proportion, fill = is_dmg)) +
#     geom_area(stat = "identity", position = "stack", alpha = 1, color = "black", linewidth = 0.2) +
#     facet_nested(~core, scales = "free")+
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#     coord_flip()+
#     theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 5/1)+
#     labs(x = "Age (y bp)", y = "Number of species", fill = "References classified as")

# png(file = "plots/microbial/weighted_number_subspecies_damage.png", width = 10, height = 10, units = "in", res = 600)
# plotA
# dev.off()


# plotA <- tax_data_dmg %>%
#   group_by(label, is_dmg) |>
#   summarise(abundance = sum(abundance), .groups = "drop") %>%
#   complete(label, is_dmg, fill = list(abundance = 0)) %>%
#   inner_join(cdata) |> 
#   mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#   mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#   ggplot(aes(x = y_bp, y = abundance, fill = is_dmg)) +
#     geom_area(stat = "identity", position = "stack", alpha = 1, color = "black", linewidth = 0.2) +
#     facet_nested(~core, scales = "free")+
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#     scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#     coord_flip()+
#     theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 5/1)+
#     labs(x = "Age (y bp)", y = "Abundance", fill = "References classified as")

# png(file = "plots/microbial/abundance_damage.png", width = 10, height = 10, units = "in", res = 600)
# plotA
# dev.off()


# plotA <- tax_data_dmg %>%
#   group_by(label, is_dmg) |>
#   summarise(abundance = sum(n_reads), .groups = "drop") %>%
#   complete(label, is_dmg, fill = list(abundance = 0)) %>%
#   inner_join(cdata) |> 
#   mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#   mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#   ggplot(aes(x = y_bp, y = abundance, fill = is_dmg)) +
#     geom_area(stat = "identity", position = "stack", alpha = 1, color = "black", linewidth = 0.2) +
#     facet_nested(~core, scales = "free")+
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#     scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#     coord_flip()+
#     theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", aspect.ratio = 5/1)+
#     labs(x = "Age (y bp)", y = "Number of reads", fill = "References classified as")

# png(file = "plots/microbial/reads_damage.png", width = 10, height = 10, units = "in", res = 600)
# plotA
# dev.off()









# png(file = "results/microbial/plots/prelim/A_b_and_fit_facet.png", width = 16, height = 8, units = "in", res = 600)
# tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = A_b, y = y_bp, color = fit))+
#         geom_point(size = 1, alpha = 0.6)+
#         geom_point(data = . %>% group_by(core, y_bp, fit) %>% summarise(mean_A_b = median(A_b), sd_A_b = sd(A_b), .groups = "drop"), 
#             aes(x = mean_A_b, fill = fit), size = 3, shape = 21, color = "black") +
    
#         # geom_errorbarh(data = . %>% group_by(core, y_bp, fit) %>% summarise(mean_A_b = mean(A_b), n = sum(!is.na(A_b)), sd_A_b = sd(A_b)/sqrt(n), .groups = "drop"), 
#         #     aes(x = mean_A_b, xmin = mean_A_b - sd_A_b, xmax = mean_A_b + sd_A_b), color = "black", height = 5, alpha = 1) +
#         # labs(x = "Damage (A_b)", y = "Age (y bp)")+
#         labs(x = "Damage (A_b)", y = "Age (y bp)", fill = "Damage fit (ccc)", color = "Damage fit (ccc)")+
#         scale_y_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#         facet_nested(~core*fit)
# dev.off()

# png(file = "results/microbial/plots/prelim/A_b_and_fit_zfit.png", width = 16, height = 8, units = "in", res = 600)
# tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     mutate(Zfit_b = ifelse(Zfit >= 2, "sig_Zfit", "non-sig_Zfit")) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = A_b, y = y_bp, color = fit))+
#         geom_point(size = 1, alpha = 0.6)+
#         # geom_point(data = . %>% group_by(core, y_bp, fit) %>% summarise(mean_A_b = mean(A_b), sd_A_b = sd(A_b), .groups = "drop"), 
#         #     aes(x = mean_A_b, fill = fit), size = 3, shape = 21, color = "black") +
    
#         # geom_errorbarh(data = . %>% group_by(core, y_bp, fit) %>% summarise(mean_A_b = mean(A_b), n = sum(!is.na(A_b)), sd_A_b = sd(A_b)/sqrt(n), .groups = "drop"), 
#         #     aes(x = mean_A_b, xmin = mean_A_b - sd_A_b, xmax = mean_A_b + sd_A_b), color = "black", height = 5, alpha = 1) +
#         # labs(x = "Damage (A_b)", y = "Age (y bp)")+
#         labs(x = "Damage (A_b)", y = "Age (y bp)", fill = "Damage fit (ccc)", color = "Damage fit (ccc)")+
#         scale_y_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#         facet_nested(~core*Zfit_b)
# dev.off()



# ## Smiley plots
# lab <- cdata %>%
#     filter(
#         core == "GeoB25202_R1",
#         # y_bp > 4.5e5 & y_bp < 4.6e5) %>%
#         y_bp == max(y_bp)) %>%
#     pull(label) 
# # lab <- lab[3]

# dat_filt <- dat1ssp %>%
#     left_join(cdata) %>%
#     filter(fit == "good") %>%
#     filter(label %in% lab) %>%
#     filter(A_b > 0.2)


# select <- dat_filt %>% filter(label == lab) %>% arrange(A_b) %>% slice(1:25) %>% pull(name) %>% unique()

# dat_filt <- dat_filt %>% filter(name %in% select)

# plots100 <- purrr::map(.x = lab, dat = dat_filt, .f = function(x, dat, ncol = 5, nrow = 2, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
#     data <- dat |>
#             filter(label == x) %>%
#             arrange(name)
#     grid_size <- calculate_plot_grid(length(data$tax_name))
#     l <- lapply(data$name, function(X) {
#     df1 <- data |>
#                 filter(name == X)
#     p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
#     p <- p + ggtitle(X)
#     return(p)
#         })
#     plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
#     }, , .progress = TRUE)

# plots100


# dat1ssp %>%
#     filter(
#         name %in% select,
#         label == lab) %>%
#     select(name, A_b, rho_c, C_b, rho_c_perm_pval, q_CI_h, c_CI_l) %>%
#     arrange(desc(A_b)) %>%
#     print(n = nrow(.))
    
#     #     mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
# #     mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))






# ### Create phyloseq object ###
# # cdata <- read.table(file = "associated_data/molecular_metadata_10022025.txt", sep = "\t", header = T)
# # tax_data_dmg <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")

# # tax_data_dmg %>%
# #     select(label, phylum, class, order, family, genus, species, subspecies, is_dmg, A_b) %>%
# #     distinct() %>%
# #     write_tsv(., file = "./results/microbial/taxonomy/dmg_ssp_in_labels.tsv.gz")

# # All
# tax_data <- tax_data_dmg

# tax_data_filt <- tax_data |>
#     select(label, subspecies, abundance)

# # let's create a phyloseq object
# tax_data_df <- tax_data_filt |>
#     pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) |>
#     as.data.frame() |>
#     column_to_rownames("subspecies")

# taxonomy_df <- tax_data |>
#     ungroup() |> 
#     select(domain, phylum, class, order, family, genus, species, subspecies) |>
#     mutate(name = subspecies) |>
#     distinct() |>
#     as.data.frame() |>
#     column_to_rownames("name")

# cdata_df <- cdata |>
#     filter(label %in% colnames(tax_data_df)) |>
#     mutate(name = label) |>
#     distinct() |>
#     as.data.frame() |>
#     column_to_rownames("name")

# # phyloseq object
# tax_data_phyloseq <- phyloseq::phyloseq(
#     otu_table(tax_data_df, taxa_are_rows = TRUE),
#     tax_table(as.matrix(taxonomy_df)),
#     sample_data(cdata_df)
# )

# # Let's save the phyloseq object in case we need it later
# saveRDS(tax_data_phyloseq, "./results/microbial/taxonomy/ps_unfiltered_all.rds")

# ps_filt <- get_st(tax_data_phyloseq, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 1e-04, norm = FALSE, glom = NULL)
# saveRDS(ps_filt, "./results/microbial/taxonomy/ps_sites_prop_all.rds")

# ps_filt <- get_st(tax_data_phyloseq, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 1e-04, norm = "median", glom = NULL)
# saveRDS(ps_filt, "./results/microbial/taxonomy/ps_stand_sites_prop_all.rds")






# ### Create phyloseq object ###
# tax_data <- tax_data_dmg

# tax_data <- tax_data %>%
#     filter(is_dmg == "Damaged")

# tax_data_filt <- tax_data |>
#     select(label, subspecies, abundance)

# # let's create a phyloseq object
# tax_data_df <- tax_data_filt |>
#     pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) |>
#     as.data.frame() |>
#     column_to_rownames("subspecies")

# taxonomy_df <- tax_data |>
#     ungroup() |>
#     select(domain, phylum, class, order, family, genus, species, subspecies) |>
#     mutate(name = subspecies) |>
#     distinct() |>
#     as.data.frame() |>
#     column_to_rownames("name")

# cdata_df <- cdata |>
#     filter(label %in% colnames(tax_data_df)) |>
#     mutate(name = label) |>
#     distinct() |>
#     as.data.frame() |>
#     column_to_rownames("name")

# # phyloseq object
# tax_data_phyloseq <- phyloseq::phyloseq(
#     otu_table(tax_data_df, taxa_are_rows = TRUE),
#     tax_table(as.matrix(taxonomy_df)),
#     sample_data(cdata_df)
# )

# # Let's save the phyloseq object in case we need it later
# saveRDS(tax_data_phyloseq, "./results/microbial/taxonomy/ps_unfiltered_dmg_ssp.rds")

# sample_sums(tax_data_phyloseq) %>%
#     density() %>% plot()

# stats_derep <- read.table(file = "taxonomic-profiling/microbial/stats/all.stats-derep-summary.tsv.gz", sep = "\t", header = T)

# stats_derep %>%
#     select(label, num_seqs) %>%
#     arrange(desc(num_seqs)) %>%
# 	left_join(cdata %>% select(label, core, y_bp)) %>%
#     filter(!is.na(core)) %>%
# 	print()

# tax_data_phyloseq_filt <- prune_samples(stats_derep %>% filter(num_seqs > 1000000) %>% pull(label), tax_data_phyloseq)
# tax_data_phyloseq_filt <- prune_taxa(taxa_sums(tax_data_phyloseq_filt) > 0, tax_data_phyloseq_filt)



# ps_filt <- get_st(tax_data_phyloseq_filt, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 0, norm = FALSE, glom = NULL)
# saveRDS(ps_filt, "./results/microbial/taxonomy/ps_sites_dmg_ssp.rds")

# ps_filt <- get_st(tax_data_phyloseq_filt, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 0, norm = "median", glom = NULL)
# saveRDS(ps_filt, "./results/microbial/taxonomy/ps_stand_sites_dmg_ssp.rds")

# ps_filt <- get_st(tax_data_phyloseq_filt, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 1e-04, norm = FALSE, glom = NULL)
# saveRDS(ps_filt, "./results/microbial/taxonomy/ps_sites_prop_dmg_ssp.rds")

# ps_filt <- get_st(tax_data_phyloseq_filt, nsites = 0.01, vcoeff = 0, mean_prop_thresh = 1e-04, norm = "median", glom = NULL)
# saveRDS(ps_filt, "./results/microbial/taxonomy/ps_stand_sites_prop_dmg_ssp.rds")












# ### Some plots ###

# tax_data_dmg <- read_tsv("results/microbial/taxonomy/dmg-summary-ssp.tsv.gz")

# # tax_data_dmg <- readRDS("./results/microbial/taxonomy/ps_unfiltered_all.rds") %>%
# #     speedyseq::psmelt() %>%
# #     as_tibble() %>%
# #     rename(abundance = Abundance) %>%
# #     separate(subspecies, into = c("subspecies", "is_dmg"), sep = ";")


# plotA <- tax_data_dmg %>%
#     filter(is_dmg == "Damaged") %>%
#     group_by(label) %>%
#     summarise(count = n_distinct(species), .groups = "drop") %>%
#     left_join(cdata) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     ggplot(., aes(x = y_bp, y = count))+
#         geom_point()+
#         geom_line()+
#         facet_nested(.~core, scales = "free")+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         coord_flip()+
#         labs(x = "Age (y b p)", y = "Number of species")+
#         theme(aspect.ratio = 5/1)

# png(file = "plots/microbial/n_sp_damage.png", width = 11, height = 10, units = "in", res = 600)
# plotA
# dev.off()


# nsp <- tax_data_dmg %>%
#     select(label, domain, phylum, class, order, family, genus, species, is_dmg, y_bp, core) %>%
#     distinct() %>%
#     group_by(label, y_bp, core, is_dmg) %>%
#     count() %>%
#     ungroup()

# nsp %>%
#     group_by(core, is_dmg) %>%
#     summarise(median_sp = round(median(n, na.rm = TRUE), 0), .groups = "drop") %>%
#     pivot_wider(names_from = is_dmg, values_from = median_sp) %>%
#     arrange(desc(core)) %>%
#     knitr::kable()


# nsp %>% 
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = y_bp, y = n, color = ifelse(n < 10, "FALSE", "TRUE")))+
#         geom_point()+
#         facet_nested(.~core*is_dmg, scales = "fixed")+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#         # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
#         coord_flip()

# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     group_by(label, core, y_bp, is_dmg, domain) %>%
#     summarise(subspecies = n(), abundance = sum(abundance), .groups = "drop") %>%
#     group_by(label, is_dmg) %>%
#     mutate(relative_abundance = abundance/sum(abundance)) %>%
#     ungroup() %>%
#     left_join(nsp) %>%
#     mutate(proportion = relative_abundance*n) %>%
#     select(label, is_dmg, domain, proportion) %>%
#     tidyr::complete(.,
#         label = unique(label),
#         is_dmg = unique(is_dmg),
#         domain = unique(domain),
#         fill = list(proportion = 0)) %>%
#     left_join(cdata %>% select(label, y_bp, core)) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = y_bp, y = proportion, fill = domain))+
#         geom_area(aes(fill = domain), stat = "identity", position = "stack", alpha = 1)+
#         geom_bar(aes(color = domain, fill = domain), stat = "identity", position = "stack")+
#         scale_fill_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73"))+
#         scale_color_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73"))+
#         facet_nested(.~core*is_dmg, scales = "free")+
#         guides(color = FALSE)+
#         # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 1e4))+
#         scale_y_continuous(
#             labels = scales::label_number(scale_cut = scales::cut_short_scale()),
#             breaks = scales::breaks_pretty(n = 3)
#         )+
#         labs(fill = "Domain", x = "Age (y b p)", y = "Numer of observed species", title = "")+
#         theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")+
#         coord_flip()

# png(file = "results/microbial/plots/prelim/n_sp_wp_domain_damage.png", width = 14, height = 10, units = "in", res = 600)
# plotA
# dev.off()



# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     group_by(label, core, y_bp, is_dmg, domain) %>%
#     summarise(subspecies = n(), abundance = sum(abundance), reads = sum(n_reads), .groups = "drop") %>%
#     select(label, is_dmg, domain, reads) %>%
#     tidyr::complete(.,
#         label = unique(label),
#         is_dmg = unique(is_dmg),
#         domain = unique(domain),
#         fill = list(reads = 0)) %>%
#     left_join(cdata) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = y_bp, y = reads))+
#         geom_area(aes(fill = domain), stat = "identity", position = "stack", alpha = 0.9)+
#         geom_bar(aes(color = domain, fill = domain), stat = "identity", position = "stack")+
#         scale_fill_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73"))+
#         scale_color_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73"))+
#         facet_nested(.~core*is_dmg, scales = "free")+
#         guides(color = FALSE)+
#         # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#         labs(fill = "Domain", x = "Age (y b p)", y = "Numer of reads", title = "Number of reads for each domain.")+
#         coord_flip()+
#         theme(axis.text.x = element_text(angle = 45, hjust = 1))

# png(file = "results/microbial/plots/prelim/n_reads_domain_damage.png", width = 15, height = 8, units = "in", res = 600)
# plotA
# dev.off()


# plotA <- tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     group_by(label, core, y_bp, is_dmg, domain) %>%
#     summarise(subspecies = n(), abundance = sum(abundance), reads = sum(n_reads), .groups = "drop") %>%
#     select(label, is_dmg, domain, abundance) %>%
#     tidyr::complete(.,
#         label = unique(label),
#         is_dmg = unique(is_dmg),
#         domain = unique(domain),
#         fill = list(abundance = 0)) %>%
#     left_join(cdata) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     ggplot(., aes(x = y_bp, y = abundance, fill = domain))+
#         geom_area(aes(fill = domain), stat = "identity", position = "stack", alpha = 0.9)+
#         geom_bar(aes(color = domain, fill = domain), stat = "identity", position = "stack")+
#         scale_fill_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73"))+
#         scale_color_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73"))+
#         facet_nested(.~core*is_dmg, scales = "free")+
#         guides(color = FALSE)+
#         # scale_x_continuous(trans = sqrt_rev_trans(), labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1000, 10000, 30000, 100000, 200000, 300000, 400000, 500000))+
#         scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 5e4))+
#         scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
#         labs(fill = "Domain", x = "Age (y b p)", y = "Abundance", title = "Abundance for each domain.")+
#         coord_flip()+
#         theme(axis.text.x = element_text(angle = 45, hjust = 1))

# png(file = "results/microbial/plots/prelim/n_abundance_domain_damage.png", width = 15, height = 8, units = "in", res = 600)
# plotA
# dev.off()


# tax_data_dmg %>%
#     filter(!grepl("GeoB25206", core)) %>%
#     group_by(label, core, y_bp, is_dmg, domain) %>%
#     summarise(subspecies = n(), abundance = sum(abundance), reads = sum(n_reads), .groups = "drop") %>%
#     select(label, is_dmg, domain, reads) %>%
#     tidyr::complete(.,
#         label = unique(label),
#         is_dmg = unique(is_dmg),
#         domain = unique(domain),
#         fill = list(reads = 0)) %>%
#     left_join(cdata) %>%
#     mutate(core = fct_relevel(core, "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
#     mutate(is_dmg = fct_relevel(is_dmg, "Non-damaged", "Damaged")) %>%
#     group_by(core, is_dmg) %>%
#     summarise(mean_nreads = round(median(reads), 0), .groups = "drop") %>%
#     pivot_wider(names_from = is_dmg, values_from = mean_nreads) %>%
#     arrange(core) %>%
#     knitr::kable()
