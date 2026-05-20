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
library(phylofactor)

setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")
source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")

cdata <- read_tsv(file = "data/metadata_v5.tsv")
tax_data <- read_tsv("results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz")
taxonomy <- tax_data %>% select(domain, phylum, class, order, family, genus, species) %>% distinct()
xrf <- read_csv("/projects/caeg/people/ngm902/apps/repos/rocs/data/combined_xrf_geochemistry_full.csv")

# Agrupar en intervalos de 0.5 cm (bin inferior). 
# Si prefieres redondear al 0.5 más cercano usa round(depth_in_core_cm*2)/2
xrf <- xrf %>%
  mutate(
    ba_ti = ba/ti,
    ca_ti = ca/ti,
    p_ti  = p/ti
  ) %>%
  select(core, depth_in_core_cm, ba_ti, ca_ti, p_ti) %>%
  mutate(depth_bin_0.5 = floor(depth_in_core_cm * 2) / 2) %>%
  select(-depth_in_core_cm) %>%
  group_by(core, depth_bin_0.5) %>%
  mutate(
    ba_ti = mean(ba_ti),
    ca_ti = mean(ca_ti),
    p_ti  = mean(p_ti),
  ) %>%
  ungroup() %>%
  distinct()


table <- tax_data %>%
  filter(
    is_dmg == "Damaged",
    core %in% c("ST8", "ST13")
    ) %>%
  left_join(xrf, by = c("core"="core", "depth_in_core_cm" = "depth_bin_0.5"))


# ########################################
# #### Phylofactorization (y_bp + core)
# ########################################

source("/projects/caeg/people/ngm902/scripts/phylofactorization-functions.r")

tree_file_bac <- "/projects/caeg/people/ngm902/apps/gtdb/207/auxillary_files/bac120_r207.sp_labels.tree"
tree_file_arc <- "/projects/caeg/people/ngm902/apps/gtdb/207/auxillary_files/ar53_r207.sp_labels.tree"

DOMAINS_KEEP <- c("d__Bacteria", "d__Archaea")
nfactors <- 9
pseudocount <- 0.5
choice_mode <- "var"
ncores <- 32

# # ----------------------------
# # 1) Build Data (species x samples) + X (predictor), using BOTH domains
# # ----------------------------
# tab_sp <- table %>%
#   filter(core == "GeoB25202_R1") %>%
#   filter(domain %in% DOMAINS_KEEP) %>%
#   filter(!is.na(species)) %>%
#   select(label, y_bp, species, abundance)

# comm <- tab_sp %>%
#   group_by(label, y_bp, species) %>%
#   summarise(abundance = sum(abundance), .groups = "drop") %>%
#   pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
#   arrange(y_bp)

# X <- comm %>% transmute(y_bp)

# mat_samp_sp <- comm %>%
#   select(-y_bp) %>%
#   column_to_rownames("label") %>%
#   as.matrix()

# mat_samp_sp <- mat_samp_sp + pseudocount
# mat_samp_sp <- mat_samp_sp / rowSums(mat_samp_sp)
# mat_samp_sp[is.na(mat_samp_sp)] <- 0

# Data <- t(mat_samp_sp)  # taxa x samples

# # ----------------------------
# # 2) Load both trees + clean tips
# # ----------------------------
# tree_bac <- read.tree(tree_file_bac)
# tree_arc <- read.tree(tree_file_arc)

# tree_bac$tip.label <- strip_outer_quotes(tree_bac$tip.label)
# tree_arc$tip.label <- strip_outer_quotes(tree_arc$tip.label)

# # sanity: ensure no duplicated tip labels across trees
# dup_tips <- intersect(tree_bac$tip.label, tree_arc$tip.label)
# if (length(dup_tips) > 0) {
#   warning("Found duplicated tip labels across bac/arc trees. Example: ", dup_tips[1])
# }

# # ----------------------------
# # 3) Combine trees into one super-tree
# # ----------------------------
# # Make a tiny 2-tip backbone and graft both trees, then drop the dummy tips.
# root <- read.tree(text = "((A:1):1,(B:1):1);")
# root$tip.label <- c("A", "B")

# tree_both <- bind.tree(root, tree_bac, where = which(root$tip.label == "A"))
# tree_both <- bind.tree(tree_both, tree_arc, where = which(tree_both$tip.label == "B"))
# tree_both <- drop.tip(tree_both, intersect(tree_both$tip.label, c("A", "B")))

# # Optional: ignore branch lengths (often safer because bac/ar not on same scale)
# tree_both$edge.length <- NULL

# # ----------------------------
# # 4) Align Data and tree
# # ----------------------------
# keep <- intersect(rownames(Data), tree_both$tip.label)
# message("Matched species (bac+arc): ", length(keep), "/", nrow(Data),
#         " (", round(100 * length(keep) / nrow(Data), 1), "%)")

# tree_both <- drop.tip(tree_both, setdiff(tree_both$tip.label, keep))
# Data <- Data[tree_both$tip.label, , drop = FALSE]

# # ----------------------------
# # 5) Run PhyloFactor (joint)
# # ----------------------------
# pf.gam <- PhyloFactor(
#   Data, tree_both, X,
#   frmla = Data ~ s(y_bp),
#   method = "gam",
#   nfactors = nfactors,
#   ncores = ncores,
#   choice = choice_mode,
#   stop.early = TRUE
# )

# pf.gam

# summary_bac <- pf_factor_summary(pf.gam)
# summary_bac


# ## Rewrite the function to work with y_bp
# plot_factor_balance <- function(PF, factor_i, flip = TRUE) {
#   m <- PF$models[[factor_i]]
#   mf <- model.frame(m)

#   depth <- mf$y_bp
#   balance <- as.numeric(model.response(mf))

#   pred <- predict(m, newdata = mf, se.fit = TRUE)
#   df <- tibble(
#     depth = depth,
#     balance = balance,
#     fit = as.numeric(pred$fit),
#     lo = as.numeric(pred$fit - 1.96 * pred$se.fit),
#     hi = as.numeric(pred$fit + 1.96 * pred$se.fit)
#   ) %>% arrange(depth)

#   p <- ggplot(df, aes(depth, balance)) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#     geom_point(alpha = 0.6) +
#     geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
#     geom_line(aes(y = fit), linewidth = 1) +
#     # scale_x_reverse() +
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 5e06, 2e4))+
#     labs(
#       x = "Depth (cm)",
#       y = "ILR balance (G1 vs G2)",
#       title = paste0("Factor ", factor_i),
#       subtitle = "Balance > 0 => G1 higher; < 0 => G2 higher"
#     )

#   if (flip) p <- p + coord_flip()
#   p
# }


# plot_all_factors(pf.gam, flip = TRUE, ncol = 3)


# B.obs <- pf.BINprojection(pf.gam, rel.abund = TRUE, prediction = FALSE, factor = factor_to_plot)

# pp <- pf.tree(
#   pf.gam,
#   bg.alpha = 1,
#   GroupList = bins(pf.gam$basis),
#   branch.length = "none",
#   top.layer = TRUE,
#   top.alpha = 0.4,
#   color.fcn = viridis::turbo
# )

# plotA_joint <- pp$ggplot

# # --- colores + orden numérico de bins
# bin_labels <- rownames(B.obs$Data)  # suelen ser "Bin 1", "Bin 2", ...
# bin_order <- bin_labels[order(as.integer(str_extract(bin_labels, "\\d+")))]

# # usa el mapping de colores que viene del árbol (si está bien):
# cols <- pp$legend$colors
# # asegúrate de que tenga nombres; si no, asigna con los bins de B.obs:
# if (is.null(names(cols))) names(cols) <- rownames(B.obs$Data)
# # reordena:
# cols <- cols[bin_order]



# y_bp <- pf.gam$X$y_bp
# ix_ybp <- order(y_bp)
# y_bp_sorted <- sort(y_bp)

# plotB_joint <- t(B.obs$Data[, ix_ybp]) %>%
#   cbind(y_bp = y_bp_sorted) %>%
#   as_tibble() %>%
#   pivot_longer(cols = -y_bp, names_to = "bin", values_to = "abundance") %>% 
#   complete(y_bp, bin, fill = list(abundance = 0)) %>%
#   mutate(
#     bin_num = as.integer(str_extract(bin, "\\d+")),
#     bin = fct_reorder(bin, bin_num)
#   ) %>%
#   select(-bin_num) %>%
#   mutate(bin = fct_relevel(bin, "Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5", "Bin 6", "Bin 7", "Bin 8", "Bin 9", "Bin 10")) %>%
#   ggplot(aes(x = y_bp, y = abundance, fill = bin)) +
#     geom_area(position = "stack", linewidth = 0.2, color = "black") +
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
#     scale_fill_manual(values = cols) +
#     labs(x = "Depth (cm)", y = "Relative abundance", fill = "Bin") +
#     theme_bw() +
#     theme(legend.position = "left") +
#     coord_flip()


# png(file = "/projects/caeg/people/ngm902/apps/repos/rocs/plots/phylofactorization_joint_bin_abundance.png", width = 20, height = 12, units = "in", res = 600)
# egg::ggarrange(plotB_joint, plotA_joint, nrow = 1, widths = c(1, 1.5), labels = c("A", "B"))
# dev.off()




# plotB_joint_facet <- t(B.obs$Data[, ix_ybp]) %>%
#   cbind(y_bp = y_bp_sorted) %>%
#   as_tibble() %>%
#   pivot_longer(cols = -y_bp, names_to = "bin", values_to = "abundance") %>% 
#   complete(y_bp, bin, fill = list(abundance = 0)) %>%
#   mutate(
#     bin_num = as.integer(str_extract(bin, "\\d+")),
#     bin = fct_reorder(bin, bin_num)
#   ) %>%
#   select(-bin_num) %>%
#   mutate(bin = fct_relevel(bin, "Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5", "Bin 6", "Bin 7", "Bin 8", "Bin 9", "Bin 10")) %>%
#   ggplot(aes(x = y_bp, y = abundance, fill = bin)) +
#     geom_area(position = "stack", linewidth = 0.2, color = "black") +
#     scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
#     scale_fill_manual(values = cols) +
#     labs(x = "Depth (cm)", y = "Relative abundance", fill = "Bin") +
#     theme_bw() +
#     facet_nested(~bin, scales = "free_x") +
#     theme(legend.position = "left") +
#     coord_flip()

# plotB_joint_facet


# # bins para los primeros K factores (consistente con B.obs cuando factor=K)
# bin_list <- bins(pf.gam$basis[, seq_len(factor_to_plot), drop = FALSE])
# # convertir a tabla: bin -> species
# bin_species <- imap_dfr(bin_list, function(idx, bin_id) {
#   tibble(
#     bin = paste0("Bin", bin_id),
#     species = pf.gam$tree$tip.label[idx]
#   )
# })

# bin_species %>% count(bin)

# bin_taxa <- bin_taxa_summary(table, bin_species,
#                             ranks = c("subspecies", "species", "genus","family","order","class", "phylum", "domain"),
#                             domain_keep = c("d__Bacteria", "d__Archaea"))



# bin_taxa %>%
#   filter(bin == "Bin9", rank %in% c("order","class", "phylum", "domain")) %>%
#   arrange(desc(prop))

# bin_taxa %>%
#   filter(bin == "Bin10", rank %in% c("order","class", "phylum", "domain")) %>%
#   arrange(desc(prop))

# bin_taxa %>%
#   filter(bin == "Bin6", rank %in% c("order","class", "phylum", "domain")) %>%
#   arrange(desc(prop))

# bin_taxa %>%
#   filter(bin == "Bin8", rank %in% c("order","class", "phylum", "domain")) %>%
#   arrange(desc(prop))


























#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################


# ----------------------------
# 1) Build Data (species x samples) + X (y_bp + core), usando TODOS los cores
# ----------------------------
tab_sp <- table %>%
  filter(domain %in% DOMAINS_KEEP) %>%
  filter(!is.na(species)) %>%
  select(label, core, ba_ti, ca_ti, p_ti, species, abundance) %>%
  mutate(
    core = factor(core),
    sample_id = paste(core, label, sep = "__")  # evita colisiones si label se repite
  )

comm <- tab_sp %>%
  group_by(sample_id, core, ba_ti, ca_ti, p_ti, species) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

# X debe tener una fila por muestra, en el MISMO orden que las filas de mat_samp_sp
X <- comm %>%
  distinct(sample_id, ba_ti, ca_ti, p_ti, core) %>%
  arrange(sample_id)

mat_samp_sp <- comm %>%
  select(-core, -ba_ti, -ca_ti, -p_ti,) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# Alinear X con las filas de la matriz (muy importante)
X <- X %>%
  filter(sample_id %in% rownames(mat_samp_sp)) %>%
  slice(match(rownames(mat_samp_sp), sample_id)) %>%
  select(ba_ti, ca_ti, p_ti, core)

# Normalización (como ya hacías)
mat_samp_sp <- mat_samp_sp + pseudocount
mat_samp_sp <- mat_samp_sp / rowSums(mat_samp_sp)
mat_samp_sp[is.na(mat_samp_sp)] <- 0

Data <- t(mat_samp_sp)  # taxa x samples

# ----------------------------
# 2) Load both trees + clean tips
# ----------------------------
tree_bac <- read.tree(tree_file_bac)
tree_arc <- read.tree(tree_file_arc)

tree_bac$tip.label <- strip_outer_quotes(tree_bac$tip.label)
tree_arc$tip.label <- strip_outer_quotes(tree_arc$tip.label)

# sanity: ensure no duplicated tip labels across trees
dup_tips <- intersect(tree_bac$tip.label, tree_arc$tip.label)
if (length(dup_tips) > 0) {
  warning("Found duplicated tip labels across bac/arc trees. Example: ", dup_tips[1])
}

# ----------------------------
# 3) Combine trees into one super-tree
# ----------------------------
# Make a tiny 2-tip backbone and graft both trees, then drop the dummy tips.
root <- read.tree(text = "((A:1):1,(B:1):1);")
root$tip.label <- c("A", "B")

tree_both <- bind.tree(root, tree_bac, where = which(root$tip.label == "A"))
tree_both <- bind.tree(tree_both, tree_arc, where = which(tree_both$tip.label == "B"))
tree_both <- drop.tip(tree_both, intersect(tree_both$tip.label, c("A", "B")))

# Optional: ignore branch lengths (often safer because bac/ar not on same scale)
tree_both$edge.length <- NULL

# ----------------------------
# 4) Align Data and tree
# ----------------------------
keep <- intersect(rownames(Data), tree_both$tip.label)
message("Matched species (bac+arc): ", length(keep), "/", nrow(Data),
        " (", round(100 * length(keep) / nrow(Data), 1), "%)")

tree_both <- drop.tip(tree_both, setdiff(tree_both$tip.label, keep))
Data <- Data[tree_both$tip.label, , drop = FALSE]





# ----------------------------
# 5) Run PhyloFactor (joint)
# ----------------------------

# frmla_use <- Data ~ s(y_bp, by = core) + core
frmla_use <- Data ~ 
  s(ba_ti, by = core) +
  s(ca_ti, by = core) +
  s(p_ti,  by = core) +
  core


# pf.gam <- PhyloFactor(
#   Data, tree_both, X,
#   frmla = frmla_use,
#   method = "gam",
#   nfactors = nfactors,
#   ncores = ncores,
#   choice = choice_mode,
#   stop.early = TRUE
# )

# save(pf.gam, file = "results/microbial/pf_gam.RData")
load("results/microbial/pf_gam.RData")

pf.gam

summary_pf <- pf_factor_summary(pf.gam)
summary_pf


############################################################
## Exploratory pipeline for PhyloFactor GAM with ratios+core
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

# ---------------------------
# User settings
# ---------------------------
RATIOS <- c("ba_ti", "ca_ti", "p_ti")
ALPHA <- 0.05
TOP_N_FACTORS <- 6
OUTDIR <- "plots/phylofactor_exploratory_ratios"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Sanity checks
# ---------------------------
stopifnot(exists("pf.gam"))
stopifnot(all(RATIOS %in% names(pf.gam$X)))
stopifnot("core" %in% names(pf.gam$X))

pf.gam$X$core <- factor(pf.gam$X$core)

# Utility: safe write
write_tsv_safe <- function(df, path) readr::write_tsv(df, path)
save_plot <- function(p, filename, w=10, h=7, dpi=300) {
  ggplot2::ggsave(filename = file.path(OUTDIR, filename), plot = p,
                  width = w, height = h, dpi = dpi, units = "in")
}

# ---------------------------
# 1) Factor ranking summary (your summary_pf)
# ---------------------------
summary_pf <- pf_factor_summary(pf.gam)

# Add simple ranking columns
summary_pf2 <- summary_pf %>%
  mutate(
    rank_p = rank(p_smooth, ties.method = "first"),
    sig = p_smooth < ALPHA
  ) %>%
  arrange(p_smooth)

print(summary_pf2, n = Inf)

write_tsv_safe(summary_pf2, file.path(OUTDIR, "factor_summary.tsv"))

# Choose factors to explore (sig first, then best p)
factors_to_explore <- summary_pf2 %>%
  filter(sig) %>%
  pull(factor)

if (length(factors_to_explore) == 0) {
  factors_to_explore <- summary_pf2 %>% slice_head(n = TOP_N_FACTORS) %>% pull(factor)
} else {
  factors_to_explore <- head(factors_to_explore, TOP_N_FACTORS)
}

message("Exploring factors: ", paste(factors_to_explore, collapse = ", "))

# ---------------------------
# 2) Extract GAM term significance per factor (smooths + parametric)
# ---------------------------
extract_gam_tables <- function(PF, factor_i) {
  m <- PF$models[[factor_i]]
  sm <- summary(m)

  # Parametric terms
  ptab <- as.data.frame(sm$p.table)
  ptab$term <- rownames(ptab)
  ptab <- as_tibble(ptab) %>%
    transmute(
      factor = factor_i,
      type = "parametric",
      term,
      estimate = Estimate,
      se = `Std. Error`,
      stat = `t value`,
      p = `Pr(>|t|)`
    )

  # Smooth terms
  stab <- as.data.frame(sm$s.table)
  stab$term <- rownames(stab)
  stab <- as_tibble(stab) %>%
    transmute(
      factor = factor_i,
      type = "smooth",
      term,
      edf = edf,
      ref_df = Ref.df,
      stat = F,
      p = `p-value`
    )

  list(parametric = ptab, smooth = stab)
}

gam_tables <- purrr::map(factors_to_explore, ~extract_gam_tables(pf.gam, .x))

param_tbl <- bind_rows(purrr::map(gam_tables, "parametric"))
smooth_tbl <- bind_rows(purrr::map(gam_tables, "smooth"))

# Parse smooth term -> ratio + core (works for "s(ba_ti):coreST13" etc.)
smooth_tbl2 <- smooth_tbl %>%
  mutate(
    ratio = case_when(
      str_detect(term, "ba_ti") ~ "ba_ti",
      str_detect(term, "ca_ti") ~ "ca_ti",
      str_detect(term, "p_ti")  ~ "p_ti",
      TRUE ~ NA_character_
    ),
    core = str_match(term, "core([^\\s\\)]+)")[,2] %>% as.character(),
    core = if_else(is.na(core), NA_character_, paste0("core", core)),
    sig = p < ALPHA
  )

write_tsv_safe(param_tbl,  file.path(OUTDIR, "gam_parametric_terms.tsv"))
write_tsv_safe(smooth_tbl2, file.path(OUTDIR, "gam_smooth_terms_parsed.tsv"))

# Quick overview: how many significant smooths per factor/ratio
smooth_overview <- smooth_tbl2 %>%
  filter(!is.na(ratio)) %>%
  group_by(factor, ratio) %>%
  summarise(
    n_terms = n(),
    n_sig = sum(sig, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_p)

print(smooth_overview, n = Inf)
write_tsv_safe(smooth_overview, file.path(OUTDIR, "smooth_overview.tsv"))

# ---------------------------
# 3) Partial effect plots (balance ~ ratio) by core for each factor
# ---------------------------
plot_factor_effect <- function(PF, factor_i, var, n = 200,
                               other_stat = c("median","mean"),
                               flip = TRUE) {
  other_stat <- match.arg(other_stat)
  m <- PF$models[[factor_i]]
  X <- PF$X

  stopifnot(var %in% names(X))
  stopifnot(all(RATIOS %in% names(X)))

  # Use global median/mean to fix other ratios
  agg_fun <- if (other_stat == "median") ~median(.x, na.rm = TRUE) else ~mean(.x, na.rm = TRUE)

  fixed <- X %>% summarise(across(all_of(RATIOS), agg_fun))

  grid <- tidyr::expand_grid(
    core = levels(X$core),
    tmp = seq(min(X[[var]], na.rm = TRUE),
              max(X[[var]], na.rm = TRUE),
              length.out = n)
  ) %>%
    rename(!!var := tmp) %>%
    mutate(
      ba_ti = if (var == "ba_ti") ba_ti else fixed$ba_ti,
      ca_ti = if (var == "ca_ti") ca_ti else fixed$ca_ti,
      p_ti  = if (var == "p_ti")  p_ti  else fixed$p_ti
    )

  pr <- predict(m, newdata = grid, se.fit = TRUE)

  df <- grid %>%
    mutate(
      fit = as.numeric(pr$fit),
      lo  = fit - 1.96 * as.numeric(pr$se.fit),
      hi  = fit + 1.96 * as.numeric(pr$se.fit)
    )

  p <- ggplot(df, aes(x = .data[[var]], y = fit, color = core, fill = core)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_bw() +
    labs(
      x = var,
      y = "Predicted ILR balance",
      title = paste0("Factor ", factor_i, " — partial effect: ", var),
      subtitle = paste0("Other ratios fixed at global ", other_stat)
    )

  if (flip) p <- p + coord_flip()
  p
}

# Produce and save partial effect plots
for (f in factors_to_explore) {
  for (v in RATIOS) {
    p <- plot_factor_effect(pf.gam, f, v, n = 200, other_stat = "median", flip = TRUE)
    save_plot(p, sprintf("factor_%02d_partial_%s.png", f, v), w = 10, h = 7, dpi = 300)
  }
}

# ---------------------------
# 4) Tree + bins projection + stacked area vs ratio (per core)
# ---------------------------
# Helper to compute bin colors from pf.tree output (like your original)
get_bin_colors <- function(pp_obj, Bobs) {
  bin_labels <- rownames(Bobs$Data)
  bin_order <- bin_labels[order(as.integer(str_extract(bin_labels, "\\d+")))]
  cols <- pp_obj$legend$colors
  if (is.null(names(cols))) names(cols) <- rownames(Bobs$Data)
  cols <- cols[bin_order]
  cols
}

# Build a stable sample_id mapping:
# We'll assume PF stores X in the same order as samples in Data columns.
# Data: taxa x samples, so samples are colnames(PF$Data) after alignment.
get_sample_meta <- function(PF) {
  # samples correspond to columns of PF$Data
  sample_id <- colnames(PF$Data)
  X <- PF$X %>% mutate(sample_id = sample_id)
  X
}

plot_bins_vs_ratio <- function(Bobs, meta, var, cols, flip = TRUE) {
  df <- t(Bobs$Data) %>%
    as_tibble(rownames = "sample_id") %>%
    left_join(meta, by = "sample_id") %>%
    pivot_longer(cols = starts_with("Bin"), names_to = "bin", values_to = "abundance") %>%
    mutate(
      bin_num = as.integer(str_extract(bin, "\\d+")),
      bin = fct_reorder(bin, bin_num)
    )

  p <- ggplot(df, aes(x = .data[[var]], y = abundance, fill = bin)) +
    geom_area(position = "stack", linewidth = 0.2, color = "black") +
    facet_wrap(~ core, scales = "free_x") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    labs(x = var, y = "Relative abundance", fill = "Bin",
         title = paste0("Bins vs ", var))

  if (flip) p <- p + coord_flip()
  p
}

meta <- get_sample_meta(pf.gam)

for (f in factors_to_explore) {
  factor_to_plot <- f

  B.obs <- pf.BINprojection(pf.gam, rel.abund = TRUE, prediction = FALSE, factor = factor_to_plot)

  pp <- pf.tree(
    pf.gam,
    bg.alpha = 1,
    GroupList = bins(pf.gam$basis),
    branch.length = "none",
    top.layer = TRUE,
    top.alpha = 0.4,
    color.fcn = viridis::turbo
  )

  plotA <- pp$ggplot
  save_plot(plotA, sprintf("factor_%02d_tree_bins.png", f), w = 12, h = 10, dpi = 300)

  cols <- get_bin_colors(pp, B.obs)

  for (v in RATIOS) {
    pB <- plot_bins_vs_ratio(B.obs, meta, v, cols, flip = TRUE)
    save_plot(pB, sprintf("factor_%02d_bins_vs_%s.png", f, v), w = 12, h = 8, dpi = 300)

    # optional: combine tree + bins plot
    comb <- pB + plotA + plot_layout(widths = c(1, 1.5))
    save_plot(comb, sprintf("factor_%02d_combined_%s.png", f, v), w = 20, h = 10, dpi = 300)
  }

  # ---------------------------
  # 5) Taxa summary per bin (requires your helper bin_taxa_summary)
  # ---------------------------
  # Build bin_species for first K factors consistent with your approach
  bin_list <- bins(pf.gam$basis[, seq_len(factor_to_plot), drop = FALSE])

  bin_species <- purrr::imap_dfr(bin_list, function(idx, bin_id) {
    tibble(
      bin = paste0("Bin", bin_id),
      species = pf.gam$tree$tip.label[idx]
    )
  })

  # This uses your existing table+helper; adjust ranks if needed
  if (exists("bin_taxa_summary")) {
    bin_taxa <- bin_taxa_summary(
      table, bin_species,
      ranks = c("subspecies","species","genus","family","order","class","phylum","domain"),
      domain_keep = c("d__Bacteria", "d__Archaea")
    )
    write_tsv_safe(bin_taxa, file.path(OUTDIR, sprintf("factor_%02d_bin_taxa.tsv", f)))
  }
}

# ---------------------------
# 6) Single report table: "most informative" effects per factor
# ---------------------------
key_effects <- smooth_tbl2 %>%
  filter(!is.na(ratio)) %>%
  group_by(factor) %>%
  arrange(p, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  arrange(factor, p)

write_tsv_safe(key_effects, file.path(OUTDIR, "key_effects_top10_per_factor.tsv"))

message("Done. Outputs written to: ", OUTDIR)













# ¿Cuántas muestras y niveles de core entraron realmente?
table(X$core)
summary(X[, c("ba_ti","ca_ti","p_ti")])
colSums(is.na(X[, c("ba_ti","ca_ti","p_ti")]))


factor_to_plot <- 9  # o el que quieras

B.obs <- pf.BINprojection(pf.gam, rel.abund = TRUE, prediction = FALSE, factor = factor_to_plot)

pp <- pf.tree(
  pf.gam,
  bg.alpha = 1,
  GroupList = bins(pf.gam$basis),
  branch.length = "none",
  top.layer = TRUE,
  top.alpha = 0.4,
  color.fcn = viridis::turbo
)

plotA_joint <- pp$ggplot



i <- 2  # prueba con 2,6,1
m <- pf.gam$models[[i]]
sm <- summary(m)

sm$p.table      # efectos paramétricos (incluye core)
sm$s.table      # smooth terms: aquí verás ba_ti/cs/ p_ti por core

i <- 6  # prueba con 2,6,1
m <- pf.gam$models[[i]]
sm <- summary(m)

sm$p.table      # efectos paramétricos (incluye core)
sm$s.table      # smooth terms: aquí verás ba_ti/cs/ p_ti por core

i <- 1  # prueba con 2,6,1
m <- pf.gam$models[[i]]
sm <- summary(m)

sm$p.table      # efectos paramétricos (incluye core)
sm$s.table      # smooth terms: aquí verás ba_ti/cs/ p_ti por core

plot_factor_effect <- function(PF, factor_i, var = c("ba_ti","ca_ti","p_ti"),
                               n = 150, flip = TRUE) {
  var <- match.arg(var)
  m <- PF$models[[factor_i]]
  X <- PF$X

  # Medianas globales para fijar los otros ratios
  meds <- X |>
    dplyr::summarise(dplyr::across(c(ba_ti, ca_ti, p_ti), ~median(.x, na.rm = TRUE)))

  # Grid: var recorre rango; otros ratios fijos
  grid <- tidyr::expand_grid(
    core = levels(X$core),
    tmp  = seq(min(X[[var]], na.rm = TRUE),
               max(X[[var]], na.rm = TRUE),
               length.out = n)
  ) |>
    dplyr::rename(!!var := tmp) |>
    dplyr::mutate(
      ba_ti = if (var == "ba_ti") ba_ti else meds$ba_ti,
      ca_ti = if (var == "ca_ti") ca_ti else meds$ca_ti,
      p_ti  = if (var == "p_ti")  p_ti  else meds$p_ti
    )

  pr <- predict(m, newdata = grid, se.fit = TRUE)

  df <- grid |>
    dplyr::mutate(
      fit = as.numeric(pr$fit),
      lo  = fit - 1.96 * as.numeric(pr$se.fit),
      hi  = fit + 1.96 * as.numeric(pr$se.fit)
    )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[var]], y = fit,
                                       color = core, fill = core)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = var,
      y = "Predicted ILR balance",
      title = paste0("Factor ", factor_i, " — efecto parcial de ", var),
      subtitle = "Otros ratios fijados en su mediana"
    )

  if (flip) p <- p + ggplot2::coord_flip()
  p
}
plot_factor_effect(pf.gam, 1, "ba_ti")














plot_factor_balance <- function(PF, factor_i, flip = TRUE) {
  m <- PF$models[[factor_i]]
  mf <- model.frame(m)

  df <- tibble(
    y_bp = mf$y_bp,
    core = mf$core,
    balance = as.numeric(model.response(mf))
  )

  pred <- predict(m, newdata = mf, se.fit = TRUE)
  df <- df %>%
    mutate(
      fit = as.numeric(pred$fit),
      lo  = as.numeric(pred$fit - 1.96 * pred$se.fit),
      hi  = as.numeric(pred$fit + 1.96 * pred$se.fit)
    ) %>%
    arrange(y_bp)

  p <- ggplot(df, aes(y_bp, balance, color = core)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.6) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = core), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1) +
    scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(
      x = "y_bp",
      y = "ILR balance (G1 vs G2)",
      title = paste0("Factor ", factor_i)
    ) +
    theme_bw()

  if (flip) p <- p + coord_flip()
  p
}












































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
     

     
     



