# Set up
library(tidyverse)
library(knitr)
library(readxl)
library(measurements)
library(ggh4x)
library(ggrepel)
library(readr)
library(scales)

source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")
setwd(dir = "/projects/caeg/people/ngm902/apps/repos/rocs")



# cdata <- read_tsv("data/metadata_v5.tsv")
# fqdup <- read_tsv("/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/fqdup/fqdup_damage_summary.tsv")


# fqdup %>%
# 	left_join(cdata) %>%
# 	mutate(core = factor(core, levels = c("NEGATIVE", "LibNTC", "ExrNTC", "ExrPTC", "LibPTC", "ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
# 	ggplot(., aes(x = core, y = d_max_5prime, color = damage_status))+
# 		geom_point()+
# 		coord_flip()

# fqdup %>%
# 	left_join(cdata) %>%
# 	filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
# 	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
# 	ggplot(., aes(x = y_bp, y = d_max_5prime, color = damage_status))+
# 		geom_point()+
# 		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
# 		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
# 		coord_flip()+
# 		labs(x = "Age (y bp)", y = "Temperature (°C)", title = "Estimated sea surface temperature", color = "Method")+
# 		facet_nested(~core, scales = "free")+
# 		scale_fill_manual(values = paired_genus)+
# 		guides(fill = guide_legend(ncol = 1))

# fqdup %>%
# 	left_join(cdata) %>%
# 	filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
# 	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
# 	ggplot(., aes(x = y_bp, y = bg_5prime, color = damage_status))+
# 		geom_point()+
# 		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
# 		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
# 		coord_flip()+
# 		labs(x = "Age (y bp)", y = "Temperature (°C)", title = "Estimated sea surface temperature", color = "Method")+
# 		facet_nested(~core, scales = "free")+
# 		scale_fill_manual(values = paired_genus)+
# 		guides(fill = guide_legend(ncol = 1))

# fqdup %>%
# 	left_join(cdata) %>%
# 	filter(core %in% c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")) %>%
# 	mutate(core = factor(core, levels = c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"))) %>%
# 	ggplot(., aes(x = y_bp, y = bg_3prime, color = damage_status))+
# 		geom_point()+
# 		scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
# 		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
# 		coord_flip()+
# 		labs(x = "Age (y bp)", y = "Temperature (°C)", title = "Estimated sea surface temperature", color = "Method")+
# 		facet_nested(~core, scales = "free")+
# 		scale_fill_manual(values = paired_genus)+
# 		guides(fill = guide_legend(ncol = 1))





# input tables
fqdup <- read_tsv("taxonomic-profiling/fqdup/fqdup_damage_summary.tsv", show_col_types = FALSE)
cdata <- read_tsv("data/metadata_v5.tsv")

core_levels <- c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")

fqdup_plot <- fqdup %>%
  left_join(cdata, by = "label") %>%
  filter(core %in% core_levels) %>%
  mutate(core = factor(core, levels = core_levels)) %>%
  mutate(
    damage_call = case_when(
      artifact ~ "Artifact",
      validated ~ "Validated",
      TRUE ~ "Indeterminate"
    ),
    library_type = factor(library_type, levels = c("double-stranded", "single-stranded", "unknown")),
    damage_status = factor(damage_status)
  )

# 1) main damage amplitude
p_dmax <- fqdup_plot %>%
  ggplot(aes(x = y_bp, y = d_max_combined, color = damage_call)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
  coord_flip() +
  facet_nested(~ core, scales = "free") +
  labs(
    x = "Age (y bp)",
    y = "Combined deamination signal (d_max_combined)",
    color = "Damage call"
  )+
  guides(color = guide_legend(nrow = 1))+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

# 2) strand asymmetry by end
fqdup_long_ends <- fqdup_plot %>%
  select(label, core, y_bp, d_max_5prime, d_max_3prime) %>%
  pivot_longer(
    cols = c(d_max_5prime, d_max_3prime),
    names_to = "end",
    values_to = "dmax"
  ) %>%
  mutate(
    end = recode(end,
      d_max_5prime = "5 prime",
      d_max_3prime = "3 prime"
    )
  )

p_ends <- fqdup_long_ends %>%
  ggplot(aes(x = y_bp, y = dmax, color = end)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
  coord_flip() +
  facet_nested(~ core, scales = "free") +
  labs(
    x = "Age (y bp)",
    y = "Terminal deamination amplitude",
    color = "End"
  )+
  guides(color = guide_legend(nrow = 1))+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

# 3) oxidation-related signal
p_oxog <- fqdup_plot %>%
  ggplot(aes(x = y_bp, y = s_oxog, color = channel_c_detected)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
	scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = scales::breaks_pretty(n = 3))+
  coord_flip() +
  facet_nested(~ core, scales = "free") +
  labs(
    x = "Age (y bp)",
    y = "Oxidation signal (s_oxog)",
    color = "Channel C detected"
  )+
  guides(color = guide_legend(nrow = 1))+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

# 4) interior damage / clustering
p_interior <- fqdup_plot %>%
  ggplot(aes(x = y_bp, y = short_log2oe, color = short_z)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
  coord_flip() +
  facet_nested(~ core, scales = "free") +
  labs(
    x = "Age (y bp)",
    y = "Interior C>T clustering (log2 O/E)",
    color = "short_z"
  ) +
  guides(color = guide_legend(nrow = 1))+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

# 5) depurination
fqdup_dep_long <- fqdup_plot %>%
  select(label, core, y_bp, enrichment_5prime = depurination_enrichment_5prime,
         enrichment_3prime = depurination_enrichment_3prime) %>%
  pivot_longer(
    cols = c(enrichment_5prime, enrichment_3prime),
    names_to = "end",
    values_to = "enrichment"
  )

p_dep <- fqdup_dep_long %>%
  ggplot(aes(x = y_bp, y = enrichment, color = end)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
  coord_flip() +
  facet_nested(~ core, scales = "free") +
  labs(
    x = "Age (y bp)",
    y = "Depurination enrichment",
    color = "End"
  )+
  guides(color = guide_legend(nrow = 1))+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

# 6) library-type confidence proxy
fqdup_bic_long <- fqdup_plot %>%
  select(label, core, y_bp, bic_ds, bic_ss, bic_bias) %>%
  pivot_longer(
    cols = c(bic_ds, bic_ss, bic_bias),
    names_to = "model",
    values_to = "bic"
  )

p_bic <- fqdup_bic_long %>%
  ggplot(aes(x = y_bp, y = bic, color = model)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = seq(0, 6e06, 2e4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
  coord_flip() +
  facet_nested(~ core, scales = "free") +
  labs(
    x = "Age (y bp)",
    y = "BIC score",
    color = "Model"
  ) 

# optional: save
ggsave("plots/damage/fqdup_dmax_combined_by_core.png", p_dmax, width = 13, height = 10, dpi = 300)
ggsave("plots/damage/fqdup_terminal_ends_by_core.png", p_ends, width = 13, height = 10, dpi = 300)
ggsave("plots/damage/fqdup_oxidation_by_core.png", p_oxog, width = 13, height = 10, dpi = 300)
ggsave("plots/damage/fqdup_interior_damage_by_core.png", p_interior, width = 13, height = 10, dpi = 300)
ggsave("plots/damage/fqdup_depurination_by_core.png", p_dep, width = 13, height = 10, dpi = 300)
ggsave("plots/damage/fqdup_bic_by_core.png", p_bic, width = 13, height = 10, dpi = 300)











# -----------------------------
# 2) Unir metadata y definir parámetros
# -----------------------------
dat <- fqdup %>%
  left_join(cdata, by = "label")

# Parámetros principales representados antes
vars_damage <- c(
  "d_max_combined",                 # deaminación terminal combinada
  "d_max_5prime",                   # deaminación 5'
  "d_max_3prime",                   # deaminación 3'
  "s_oxog",                         # oxidación
  "short_log2oe",                   # daño interno
  "depurination_enrichment_5prime", # depurinación 5'
  "depurination_enrichment_3prime"  # depurinación 3'
)

# Opcional: restringir a cores concretos
core_levels <- c("ST5", "ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2")

dat <- dat %>%
  filter(core %in% core_levels) %>%
  mutate(core = factor(core, levels = core_levels))

# -----------------------------
# 3) Función para correlaciones por core
# -----------------------------
cor_by_core <- function(df, vars, method = "spearman") {
  
  # Mantener solo columnas necesarias
  df2 <- df %>%
    select(core, all_of(vars))
  
  # Todas las combinaciones de variables
  pairs <- expand.grid(var1 = vars, var2 = vars, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    filter(var1 != var2)
  
  # Calcular por core
  res <- df2 %>%
    group_by(core) %>%
    group_modify(~{
      map_dfr(seq_len(nrow(pairs)), function(i) {
        v1 <- pairs$var1[i]
        v2 <- pairs$var2[i]
        
        tmp <- .x %>%
          select(all_of(c(v1, v2))) %>%
          drop_na()
        
        n <- nrow(tmp)
        
        if (n < 3) {
          tibble(
            var1 = v1,
            var2 = v2,
            n = n,
            estimate = NA_real_,
            p_value = NA_real_
          )
        } else {
          ct <- suppressWarnings(
            cor.test(tmp[[v1]], tmp[[v2]], method = method, exact = FALSE)
          )
          
          tibble(
            var1 = v1,
            var2 = v2,
            n = n,
            estimate = unname(ct$estimate),
            p_value = ct$p.value
          )
        }
      })
    }) %>%
    ungroup() %>%
    group_by(core) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      sig = case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    ungroup()
  
  res
}

cors <- cor_by_core(dat, vars_damage, method = "spearman")

# Guardar tabla larga
readr::write_tsv(cors, "taxonomic-profiling/fqdup/fqdup_correlations_by_core.tsv")

# -----------------------------
# 4) Heatmap por core
# -----------------------------
# Para no duplicar información, nos quedamos con triángulo superior
vars_order <- vars_damage

cors_plot <- cors %>%
  mutate(
    var1 = factor(var1, levels = vars_order),
    var2 = factor(var2, levels = vars_order)
  ) %>%
  filter(as.integer(var1) < as.integer(var2))

p_heat <- cors_plot %>%
  ggplot(aes(x = var1, y = var2, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(sprintf("%.2f", estimate), sig)), size = 3) +
  facet_wrap(~core) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  labs(
    x = "",
    y = "",
    fill = "Spearman rho"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave("plots/damage/fqdup_correlations_by_core_heatmap.png", p_heat, width = 13, height = 9, dpi = 300)

# # -----------------------------
# # 5) Tabla ancha opcional por core
# # -----------------------------
# # Útil para inspección manual
# cors_wide <- cors %>%
#   mutate(pair = paste(var1, var2, sep = " vs ")) %>%
#   select(core, pair, estimate, p_adj, n) %>%
#   pivot_wider(names_from = pair, values_from = c(estimate, p_adj, n))

# readr::write_tsv(cors_wide, "taxonomic-profiling/fqdup/fqdup_correlations_by_core_wide.tsv")

# # -----------------------------
# # 6) Scatterplots solo para correlaciones significativas
# # -----------------------------
# sig_pairs <- cors %>%
#   filter(!is.na(p_adj), p_adj < 0.05) %>%
#   distinct(core, var1, var2)

# plot_pair_core <- function(df, core_name, xvar, yvar) {
#   df %>%
#     filter(core == core_name) %>%
#     ggplot(aes(x = .data[[xvar]], y = .data[[yvar]])) +
#     geom_point(alpha = 0.8) +
#     geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
#     labs(
#       title = paste(core_name, "-", xvar, "vs", yvar),
#       x = xvar,
#       y = yvar
#     ) +
#     theme_bw()
# }

# # Ejemplo: guardar todos los plots significativos
# dir.create("plots/damage/fqdup_cor_plots", showWarnings = FALSE)

# pwalk(
#   sig_pairs,
#   function(core, var1, var2) {
#     p <- plot_pair_core(dat, core, var1, var2)
#     ggsave(
#       filename = file.path("plots/damage/fqdup_cor_plots", paste0(core, "__", var1, "__", var2, ".png")),
#       plot = p,
#       width = 5,
#       height = 4,
#       dpi = 300
#     )
#   }
# )
