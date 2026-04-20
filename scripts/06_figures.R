#!/usr/bin/env Rscript
# 06_figures.R — Comparison figures: new vs old results + TEA vs EMP

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(here("config.R"))

# Path to pre-computed original analysis outputs used for new-vs-old comparisons.
# Update this if the original results are moved; it is intentionally not in
# config.R because it is system-specific and read-only reference data.
OLD_WGCNA <- "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1"
FIGS      <- RESULTS$figures

theme_ms <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── Load new results ───────────────────────────────────────────────────────────

log_msg("Loading new results...")
meta      <- fread(UPSTREAM$metadata)
new_MEs   <- fread(file.path(RESULTS$stage1, "wgcna", "module_eigengenes.tsv"))
new_mods  <- fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
new_pres  <- fread(file.path(RESULTS$stage1, "wgcna", "preservation.tsv"))
new_hmm   <- fread(file.path(RESULTS$hmm, "hmm_states.tsv"))
new_bic   <- fread(file.path(RESULTS$hmm, "bic_comparison.tsv"))
new_emp   <- fread(file.path(RESULTS$emp, "emp_sap_per_sample.tsv"))
new_tea   <- fread(file.path(RESULTS$tea, "tea_indices_per_sample.tsv"))
tea_corr  <- fread(file.path(RESULTS$tea, "tea_vs_emp_correlations.tsv"))
tea_clim  <- fread(file.path(RESULTS$tea, "tea_climate_models.tsv"))

# Merge metadata into new MEs (temp_complete loaded for potential future use)
new_trait <- merge(new_MEs, meta[, .(label, core, y_bp, mis, temp_complete)],
                   by.x = "sample", by.y = "label")
new_trait[, age_kyr := y_bp / 1000]
setnames(new_trait, "mis", "d18O")

# ── Load old results ───────────────────────────────────────────────────────────

log_msg("Loading old results...")
old_MEs  <- fread(file.path(OLD_WGCNA, "wgcna_prokaryotes", "module_eigengenes.tsv"))
old_mods <- fread(file.path(OLD_WGCNA, "wgcna_prokaryotes", "module_assignments.tsv"))
old_hmm  <- fread(file.path(OLD_WGCNA, "hmm_states", "hmm_states.tsv"))
old_emp  <- fread(file.path(OLD_WGCNA, "emp_sap", "emp_sap_per_sample.tsv"))

old_trait <- merge(old_MEs, meta[, .(label, core, y_bp, mis, temp_complete)],
                   by.x = "sample", by.y = "label")
old_trait[, age_kyr := y_bp / 1000]
setnames(old_trait, "mis", "d18O")

# ── Fig 1: Module composition comparison ─────────────────────────────────────

log_msg("Fig 1: module sizes...")

new_sz <- new_mods[, .(n = .N, run = "new"), by = module][module != "grey"]
old_sz <- old_mods[, .(n = .N, run = "old"), by = module][module != "grey"]
mod_sz <- rbind(new_sz, old_sz)

p_mods <- ggplot(mod_sz, aes(x = reorder(module, -n), y = n, fill = run)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(new = "#2196F3", old = "#FF9800"), name = NULL) +
  labs(title = "A. Module sizes (new vs old)", x = "Module", y = "Taxa") +
  theme_ms

# ── Fig 2: Module preservation ───────────────────────────────────────────────

log_msg("Fig 2: preservation...")

pres_plot <- new_pres[!module %in% c("gold", "grey")]
p_pres <- ggplot(pres_plot, aes(x = reorder(module, -Zsummary), y = Zsummary, fill = preserved)) +
  geom_col() +
  geom_hline(yintercept = c(2, 10), linetype = "dashed", colour = c("orange", "red")) +
  scale_fill_manual(values = c(strong = "#4CAF50", moderate = "#FF9800", weak = "#F44336"),
                    name = "Preservation") +
  labs(title = "B. Module preservation (GeoB R1 → R2)",
       x = "Module", y = "Zsummary") +
  theme_ms

# ── Fig 3: Eigengenes vs d18O (new) ──────────────────────────────────────────

log_msg("Fig 3: eigengenes vs d18O...")

me_cols_new <- grep("^ME(?!grey)", names(new_trait), value = TRUE, perl = TRUE)
me_long <- melt(new_trait[, c("sample", "age_kyr", "d18O", "core", me_cols_new), with = FALSE],
                id.vars = c("sample", "age_kyr", "d18O", "core"),
                variable.name = "module", value.name = "eigengene")
me_long[, module := sub("^ME", "", module)]

p_me_d18o <- ggplot(me_long, aes(x = d18O, y = eigengene)) +
  geom_point(aes(colour = core), alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.7, se = TRUE) +
  facet_wrap(~ module, scales = "free_y", nrow = 2) +
  scale_colour_brewer(palette = "Set2", name = "Core") +
  labs(title = "C. Module eigengenes vs δ¹⁸O (new)", x = "δ¹⁸O (‰)", y = "Eigengene") +
  theme_ms + theme(legend.position = "top")

# ── Fig 4: HMM states BIC + state composition ────────────────────────────────

log_msg("Fig 4: HMM...")

p_bic <- ggplot(new_bic, aes(x = K, y = BIC)) +
  geom_line() + geom_point(size = 3) +
  geom_point(data = new_bic[K == 5], colour = "red", size = 4) +
  labs(title = "D. HMM BIC (K=5 selected)", x = "Number of states", y = "BIC") +
  theme_ms

# State × climate
new_hmm <- new_hmm
state_d18o <- new_hmm[, .(mean_d18O = mean(d18O, na.rm = TRUE), n = .N), by = label]

p_state_d18o <- ggplot(new_hmm, aes(x = label, y = d18O, fill = label)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "E. HMM states vs δ¹⁸O", x = "State", y = "δ¹⁸O (‰)") +
  theme_ms

# ── Fig 5: EMP new vs old correlation ────────────────────────────────────────

log_msg("Fig 5: EMP comparison...")

emp_both <- merge(
  new_emp[, .(sample, EMP_new = EMP_scaled)],
  old_emp[, .(sample, EMP_old = EMP_scaled)],
  by = "sample"
)
r_emp <- cor(emp_both$EMP_new, emp_both$EMP_old, use = "complete.obs")

p_emp_comp <- ggplot(emp_both, aes(x = EMP_old, y = EMP_new)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", colour = "steelblue") +
  labs(title = sprintf("F. EMP: new vs old (r = %.3f)", r_emp),
       x = "EMP old (scaled)", y = "EMP new (scaled)") +
  theme_ms

# ── Fig 6: EMP and TEA vs d18O ────────────────────────────────────────────────

log_msg("Fig 6: EMP and TEA vs d18O...")

full_dt <- Reduce(function(a, b) merge(a, b, by = "sample", all.x = TRUE),
                  list(new_emp[, .(sample, EMP_scaled, SAP)],
                       new_tea[, .(sample, OAP, DCI, MII, SRPI)],
                       meta[, .(sample = label, d18O = mis, core, age_kyr = y_bp/1000)]))

idx_long <- melt(full_dt[, .(sample, d18O, core, EMP = EMP_scaled, OAP, MII, DCI, SRPI)],
                 id.vars = c("sample", "d18O", "core"),
                 variable.name = "index", value.name = "value")
idx_long <- idx_long[is.finite(value)]

# z-score within index so all indices are on the same axis for visual comparison.
# This removes absolute magnitude differences (e.g. OAP is in mV, EMP in J/sample)
# but preserves direction and relative variation — do not interpret z-scores as
# equivalent effect sizes.
idx_long[, value_z := scale(value), by = index]

p_tea_d18o <- ggplot(idx_long, aes(x = d18O, y = value_z)) +
  geom_point(aes(colour = core), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = FALSE) +
  facet_wrap(~ index, nrow = 1, scales = "free_y") +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  labs(title = "G. EMP and TEA indices vs δ¹⁸O (z-scored)",
       x = "δ¹⁸O (‰)", y = "z-score") +
  theme_ms

# ── Fig 7: TEA correlation heatmap ────────────────────────────────────────────

log_msg("Fig 7: TEA correlation heatmap...")

# Build a square correlation matrix from the upper-triangle tea_corr table.
# tea_vs_emp_correlations.tsv contains only index1 < index2 pairs; we mirror
# them and add the diagonal (r=1). Both axes use unique(tea_corr$index1), so if
# an index appears only as index2 in the source file it will be absent here —
# verify the source file contains all indices on the index1 side.
corr_sq <- rbind(
  tea_corr[, .(index1, index2, r)],
  tea_corr[, .(index1 = index2, index2 = index1, r)],
  data.table(index1 = unique(tea_corr$index1),
             index2 = unique(tea_corr$index1),
             r = 1)
)
idx_order <- c("EMP", "OAP", "MII", "DCI", "SRPI", "SAP")
idx_order <- intersect(idx_order, unique(corr_sq$index1))
corr_sq[, index1 := factor(index1, levels = idx_order)]
corr_sq[, index2 := factor(index2, levels = rev(idx_order))]

p_corr <- ggplot(corr_sq[!is.na(index1) & !is.na(index2)],
                 aes(x = index1, y = index2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "H. Index pairwise correlations") +
  theme_ms + theme(axis.title = element_blank())

# ── Assemble and save ─────────────────────────────────────────────────────────

log_msg("Assembling figures...")

# Figure dimensions target A4 half-page (280 × 200 mm) or quarter-page (200 × 100 mm)
# at 150 dpi for screen review and 300 dpi for print (change dpi= for publication).
# Page 1: WGCNA
p1 <- (p_mods | p_pres) / p_me_d18o +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(title = "WGCNA Results: New vs Old", theme = theme(plot.title = element_text(size = 12)))

ggsave(file.path(FIGS, "fig1_wgcna.pdf"),   p1, width = 280, height = 200, units = "mm")
ggsave(file.path(FIGS, "fig1_wgcna.png"),   p1, width = 280, height = 200, units = "mm", dpi = 150)

# Page 2: HMM
p2 <- (p_bic | p_state_d18o) +
  plot_annotation(title = "HMM States", theme = theme(plot.title = element_text(size = 12)))

ggsave(file.path(FIGS, "fig2_hmm.pdf"),  p2, width = 200, height = 100, units = "mm")
ggsave(file.path(FIGS, "fig2_hmm.png"),  p2, width = 200, height = 100, units = "mm", dpi = 150)

# Page 3: EMP + TEA
p3 <- (p_emp_comp | p_corr) / p_tea_d18o +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "EMP / TEA Comparison", theme = theme(plot.title = element_text(size = 12)))

ggsave(file.path(FIGS, "fig3_emp_tea.pdf"), p3, width = 280, height = 200, units = "mm")
ggsave(file.path(FIGS, "fig3_emp_tea.png"), p3, width = 280, height = 200, units = "mm", dpi = 150)

log_msg("Figures saved to ", FIGS)

# ── Print summary table ────────────────────────────────────────────────────────

cat("\n=== RESULTS COMPARISON ===\n")
cat(sprintf("%-30s %10s %10s\n", "Metric", "New", "Old"))
cat(strrep("-", 52), "\n")
cat(sprintf("%-30s %10d %10d\n", "Samples (stage1)",        nrow(new_MEs), nrow(old_MEs)))
cat(sprintf("%-30s %10d %10d\n", "WGCNA taxa",              nrow(new_mods), nrow(old_mods)))
cat(sprintf("%-30s %10d %10d\n", "Non-grey modules",        new_mods[module!="grey", uniqueN(module)], old_mods[module!="grey", uniqueN(module)]))
new_emp_summary <- fread(file.path(RESULTS$emp, "emp_sap_summary.tsv"))
new_emp_r <- new_emp_summary[metric == "EMP", cor_d18O]
old_emp_m <- merge(old_emp, meta[, .(label, d18O = mis)], by.x = "sample", by.y = "label")
old_emp_r <- cor(old_emp_m$EMP, old_emp_m$d18O, use = "complete.obs")
cat(sprintf("%-30s %10.3f %10.3f\n", "EMP ~ d18O (r)", new_emp_r, old_emp_r))
cat(sprintf("%-30s %10s %10s\n", "Significant TEA (FDR<0.05)",
            paste(tea_clim[!is.na(fdr) & fdr < 0.05, index], collapse=","), "—"))
cat("\n")
