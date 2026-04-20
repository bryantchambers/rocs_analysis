#!/usr/bin/env Rscript
# =============================================================================
# 13_fig_module_composition.R - Module Composition Figures (ds4_size15)
# =============================================================================
#
# Updated for 13 modules with ds4_size15 parameters
# Uses M1-M13 labels instead of color names
#
# =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(ggnewscale)

source("/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/lib_figure_utils.R")

config <- list(
  base_dir = "/projects/caeg/people/kbd606/scratch/mateu-rocs",
  results_dir = "/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1",
  class_dir = "/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/classification",
  output_dir = "/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/figures",
  target_cores = c("GeoB25202_R1", "GeoB25202_R2", "ST13", "ST8"),
  outlier_samples = c("LV3003046968"),
  mis_file = "/projects/caeg/people/kbd606/scratch/mateu-rocs/results/network_global/mis_stage_boundaries.tsv",
  lr04_file = "/projects/fernandezguerra/people/ngm902/ROCS/associated_data/Lisiecki2005_copy.txt"
)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# 7 modules from CLR baseline (ordered by size)
module_order <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black")

# M1-M7 labels
module_labels <- setNames(paste0("M", 1:7), module_order)

# Climate-linked modules (6 significant: 4 glacial_UP, 2 glacial_DOWN)
# glacial_UP: blue, brown, turquoise, yellow
# glacial_DOWN: black, red
climate_modules <- c("blue", "brown", "turquoise", "yellow", "black", "red")

# Colors for modules
module_colors <- c(
  "turquoise" = "#1B9E77", "blue" = "#377EB8", "brown" = "#A65628",
  "yellow" = "#FFFF33", "green" = "#66A61E", "red" = "#E7298A",
  "black" = "#666666", "pink" = "#FB8072", "magenta" = "#984EA3",
  "purple" = "#7570B3", "greenyellow" = "#B3DE69", "tan" = "#D2B48C", "salmon" = "#FA8072"
)

main <- function() {
  log_msg("=", strrep("=", 60))
  log_msg("Module Composition Figures (CLR baseline)")
  log_msg("=", strrep("=", 60))

  # ---------------------------------------------------------------------------
  # Load data
  # ---------------------------------------------------------------------------
  log_msg("\nLoading data...")

  meta <- fread(file.path(config$base_dir, "data/metadata_v4.txt"))
  meta <- meta[core %in% config$target_cores & !label %in% config$outlier_samples]
  meta[, age_kyr := y_bp / 1000]
  meta <- meta[age_kyr <= 150]
  log_msg(sprintf("  Stage 1 samples: %d", nrow(meta)))

  modules <- fread(file.path(config$results_dir, "wgcna_prokaryotes/module_assignments.tsv"))
  classification <- fread(file.path(config$class_dir, "prokaryote_function_assigned.tsv"))

  mis_bounds <- fread(config$mis_file)
  setnames(mis_bounds, c("stage", "start_kyr", "end_kyr", "climate"), c("MIS", "start", "end", "climate"))
  mis_bounds <- mis_bounds[end <= 160]

  lr04 <- fread(config$lr04_file)
  setnames(lr04, c("Time_ka", "Benthic_d18O", "Standard_error"), c("age_kyr", "d18O", "se"))
  lr04_150 <- lr04[age_kyr <= 150]

  mod_class <- merge(modules, classification[, .(taxon, functional_group, signal_source)],
                     by = "taxon", all.x = TRUE)
  mod_class[is.na(functional_group), functional_group := "Unknown"]

  log_msg(sprintf("  Total modules: %d", length(unique(modules[module != "grey"]$module))))
  log_msg(sprintf("  Climate-linked: %d (%s)", length(climate_modules),
                  paste(module_labels[climate_modules], collapse=", ")))

  # ---------------------------------------------------------------------------
  # Load raw prokaryote data
  # ---------------------------------------------------------------------------
  log_msg("\nLoading raw prokaryote reads...")
  prok_raw <- fread(file.path(config$base_dir, "results/microbial/taxonomy/dmg-summary-ssp_selected.tsv.gz"))
  prok_raw <- prok_raw[label %in% meta$label & core %in% config$target_cores & is_dmg == "Damaged"]
  prok_raw <- merge(prok_raw[, .(subspecies, label, core, n_reads)],
                    mod_class[, .(taxon, module, functional_group)],
                    by.x = "subspecies", by.y = "taxon", all.x = FALSE)
  prok_raw <- merge(prok_raw, meta[, .(label, age_kyr)], by = "label")

  # ===========================================================================
  # PANEL: LR04 δ¹⁸O Reference
  # ===========================================================================
  log_msg("\nCreating LR04 reference panel...")

  mis_labels <- mis_bounds[start < 150]
  mis_labels[, label_x := (start + pmin(end, 150)) / 2]
  d18o_range <- range(lr04_150$d18O, na.rm = TRUE)

  p_climate <- ggplot() +
    geom_rect(data = mis_bounds[start < 150],
              aes(xmin = start, xmax = pmin(end, 150), ymin = -Inf, ymax = Inf, fill = climate),
              alpha = 0.6) +
    scale_fill_manual(values = mis_colors, guide = "none") +
    geom_text(data = mis_labels,
              aes(x = label_x, y = d18o_range[1] - 0.1, label = MIS),
              size = 3, color = "grey30", fontface = "bold") +
    geom_line(data = lr04_150, aes(x = age_kyr, y = d18O), color = "black", linewidth = 0.5) +
    scale_x_reverse(limits = c(150, 0), breaks = seq(150, 0, -25), expand = c(0.01, 0)) +
    scale_y_reverse() +
    labs(x = NULL, y = expression(delta^{18}*"O (‰)"), title = "LR04 δ¹⁸O (benthic stack)") +
    theme_nature(base_size = 9) +
    theme(axis.text.x = element_blank(), plot.margin = margin(5, 5, 0, 5))

  # ===========================================================================
  # PLOT A: Module abundance through time (climate-linked only)
  # ===========================================================================
  log_msg("Creating Panel A: Climate-linked module abundance...")

  module_sample <- prok_raw[module %in% climate_modules,
                             .(n_reads = sum(n_reads)),
                             by = .(module, label, core, age_kyr)]
  module_sample[, total_reads := sum(n_reads), by = label]
  module_sample[, rel_abund := n_reads / total_reads]

  module_sample[, age_bin := floor(age_kyr / 10) * 10 + 5]
  module_binned <- module_sample[, .(n_reads = sum(n_reads)), by = .(module, core, age_bin)]
  module_binned[, total_reads := sum(n_reads), by = .(core, age_bin)]
  module_binned[, rel_abund := n_reads / total_reads]

  all_combos <- CJ(module = climate_modules, core = unique(module_binned$core),
                   age_bin = unique(module_binned$age_bin))
  module_binned <- merge(all_combos, module_binned, by = c("module", "core", "age_bin"), all.x = TRUE)
  module_binned[is.na(rel_abund), rel_abund := 0]

  module_binned[, module := factor(module, levels = rev(climate_modules))]
  module_binned[, module_label := module_labels[as.character(module)]]
  module_binned[, core_label := factor(core_labels[core],
                                        levels = c("GeoB25202-R1", "GeoB25202-R2", "ST13", "ST8"))]

  p_abundance <- ggplot() +
    geom_rect(data = mis_bounds[start < 150],
              aes(xmin = start, xmax = pmin(end, 150), ymin = 0, ymax = 1, fill = climate), alpha = 0.6) +
    scale_fill_manual(values = mis_colors, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_area(data = module_binned, aes(x = age_bin, y = rel_abund, fill = module),
              position = "stack", alpha = 0.9) +
    scale_fill_manual(values = module_colors[climate_modules],
                      labels = module_labels[climate_modules],
                      guide = guide_legend(reverse = TRUE)) +
    facet_wrap(~ core_label, ncol = 1) +
    scale_x_reverse(limits = c(150, 0), breaks = seq(150, 0, -25), expand = c(0.01, 0)) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 1)) +
    labs(x = "Age (ka)", y = "Relative abundance",
         title = "A. Climate-Linked Module Abundance (6 validated modules)",
         fill = "Module") +
    theme_nature(base_size = 9) +
    theme(legend.key.size = unit(0.3, "cm"))

  # ===========================================================================
  # PLOT B: Functional composition per module
  # ===========================================================================
  log_msg("Creating Panel B: Category composition by module...")

  category_counts <- mod_class[module %in% climate_modules & !is.na(functional_group),
                                .N, by = .(module, functional_group)]
  category_counts[, prop := N / sum(N), by = module]
  category_counts[, module := factor(module, levels = climate_modules)]
  category_counts[, module_label := module_labels[as.character(module)]]

  cat_order <- category_counts[, .(total = sum(N)), by = functional_group][order(-total)]$functional_group
  category_counts[, functional_group := factor(functional_group, levels = rev(cat_order))]

  p_category <- ggplot(category_counts, aes(x = module_label, y = prop, fill = functional_group)) +
    geom_col(position = "stack", width = 0.75) +
    scale_fill_manual(values = prok_colors, guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.02)) +
    labs(x = "Module", y = "Proportion of taxa",
         title = "B. Functional Composition (Climate-Linked Modules)",
         fill = "Category") +
    theme_nature(base_size = 9) +
    theme(legend.key.size = unit(0.3, "cm"))

  # ===========================================================================
  # PLOT C: All modules summary (climate vs non-climate)
  # ===========================================================================
  log_msg("Creating Panel C: All modules overview...")

  all_mod_counts <- mod_class[module != "grey", .N, by = module]
  all_mod_counts[, module := factor(module, levels = module_order)]
  all_mod_counts[, module_label := module_labels[as.character(module)]]
  all_mod_counts[, is_climate := module %in% climate_modules]

  p_overview <- ggplot(all_mod_counts, aes(x = module_label, y = N, fill = is_climate)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = N), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
                      labels = c("TRUE" = "Climate-linked", "FALSE" = "Not validated"),
                      name = NULL) +
    labs(x = "Module", y = "Number of taxa",
         title = "C. Module Sizes (6/13 climate-linked)") +
    theme_nature(base_size = 9) +
    theme(legend.position = "top")

  # ===========================================================================
  # Combine and save
  # ===========================================================================
  log_msg("\nCombining and saving...")
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  p_combined <- p_climate / p_abundance / (p_category | p_overview) +
    plot_layout(heights = c(0.5, 1.5, 1))

  ggsave(file.path(config$output_dir, "fig_module_composition.pdf"), p_combined,
         width = 300, height = 360, units = "mm")
  ggsave(file.path(config$output_dir, "fig_module_composition.png"), p_combined,
         width = 300, height = 360, units = "mm", dpi = 300)
  log_msg("  Saved: fig_module_composition.pdf/png")

  # Export data
  fwrite(category_counts, file.path(config$output_dir, "fig_module_category_data.tsv"), sep = "\t")

  # Print summary
  log_msg("\n=== Climate-Linked Modules Summary ===")
  summary_dt <- mod_class[module %in% climate_modules, .N, by = module][order(-N)]
  summary_dt[, module_label := module_labels[module]]
  print(summary_dt[, .(module_label, taxa = N)])

  log_msg("\nDone!")
}

main()
