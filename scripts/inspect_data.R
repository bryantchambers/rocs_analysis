#!/usr/bin/env Rscript
# scripts/inspect_data.R - Write a compact, current repository data summary.
#
# The output DATA_SUMMARY.md is intended as working memory for future analysis:
# it records what each major output looks like, how samples/taxa align, and
# which generated tables should be trusted before rerunning expensive steps.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

output_file <- here("DATA_SUMMARY.md")

fmt_num <- function(x, digits = 4) {
  if (length(x) == 0 || is.na(x[1])) return("NA")
  formatC(as.numeric(x), digits = digits, format = "fg", flag = "#")
}

fmt_bytes <- function(x) {
  x <- as.numeric(x)
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  if (!any(ok)) return(out)
  units <- c("B", "KB", "MB", "GB", "TB")
  idx <- pmax(1, pmin(length(units), floor(log(pmax(x[ok], 1), 1024)) + 1))
  vals <- x[ok] / (1024 ^ (idx - 1))
  out[ok] <- paste0(formatC(vals, digits = 3, format = "fg", flag = "#"), " ", units[idx])
  out
}

exists_mark <- function(path) if (file.exists(path)) "present" else "missing"

safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  fread(path, ...)
}

md_table <- function(dt, max_rows = Inf) {
  if (is.null(dt) || nrow(dt) == 0) {
    cat("_None_\n\n")
    return(invisible(NULL))
  }
  dt <- as.data.table(dt)
  if (is.finite(max_rows)) dt <- head(dt, max_rows)
  cols <- names(dt)
  cat("|", paste(cols, collapse = " | "), "|\n", sep = "")
  cat("|", paste(rep("---", length(cols)), collapse = " | "), "|\n", sep = "")
  for (i in seq_len(nrow(dt))) {
    vals <- vapply(dt[i], function(x) {
      x <- as.character(x)
      ifelse(is.na(x) | x == "", "NA", gsub("\\|", "/", x))
    }, character(1))
    cat("|", paste(vals, collapse = " | "), "|\n", sep = "")
  }
  cat("\n")
}

schema_row <- function(path, label = path) {
  dt <- safe_fread(path, nrows = 0)
  if (is.null(dt)) {
    return(data.table(file = label, status = "missing", rows = NA_integer_,
                      cols = NA_integer_, bytes = NA_character_, columns = NA_character_))
  }
  n_rows <- nrow(fread(path, select = 1, showProgress = FALSE))
  data.table(
    file = label,
    status = "present",
    rows = n_rows,
    cols = length(names(dt)),
    bytes = fmt_bytes(file.info(path)$size),
    columns = paste(names(dt), collapse = ", ")
  )
}

results_inventory <- function() {
  if (!dir.exists(RESULTS_DIR)) return(data.table())
  rel <- list.files(RESULTS_DIR, recursive = TRUE, full.names = FALSE)
  full <- file.path(RESULTS_DIR, rel)
  info <- file.info(full)
  keep <- !is.na(info$isdir) & !info$isdir
  data.table(
    file = file.path("results", rel[keep]),
    bytes_num = as.numeric(info$size[keep])
  )[, ext := sub("^.*\\.", "", file)]
}

output_group_summary <- function(inv) {
  if (is.null(inv) || nrow(inv) == 0) return(data.table())
  grouped <- copy(inv)
  grouped[, group := fifelse(
    grepl("^results/microbial/sourcetracker/[^/]+/", file),
    sub("^(results/microbial/sourcetracker/[^/]+).*", "\\1", file),
    fifelse(
      grepl("^results/microbial/sourcetracker/", file),
      "results/microbial/sourcetracker",
      fifelse(
        grepl("^results/microbial/damage/", file),
        "results/microbial/damage",
        fifelse(
          grepl("^results/microbial/taxonomy/", file),
          "results/microbial/taxonomy",
          sub("^(results/[^/]+).*", "\\1", file)
        )
      )
    )
  )]
  out <- grouped[, .(
    files = .N,
    total_bytes_num = sum(bytes_num, na.rm = TRUE),
    largest_file = file[which.max(bytes_num)],
    largest_bytes_num = max(bytes_num, na.rm = TRUE)
  ), by = group][order(-total_bytes_num)]
  out[, total_size := fmt_bytes(total_bytes_num)]
  out[, largest_size := fmt_bytes(largest_bytes_num)]
  out[, c("total_bytes_num", "largest_bytes_num") := NULL]
  out[]
}

pipeline_steps <- data.table(
  step = c("01", "02", "03", "04", "05", "06", "06b", "07", "07b", "08",
           "09", "10", "11", "12", "13", "14", "15", "16"),
  script = c("01_data_prep.R", "02_wgcna.R", "03_hmm_states.R", "04_emp.R",
             "05_tea_vs_emp.R", "06_figures.R", "06b_bryantfigures.R",
             "07_taxon_importance.R", "07b_taxon_importance_fuzzy.R",
             "08_network_statistics.R", "09_driver_integration.R",
             "10_climate_sensitivity.R", "11_state_networks.R",
             "12_functional_linkage.R", "13_state_transition_network.R",
             "14_driver_quadrants.R", "15_state_functional_breakdown.R",
             "16_final_story_visualization.R"),
  output_role = c("Stage-1 filtering, CLR/VST object, metadata",
                  "Consensus WGCNA modules and preservation",
                  "HMM ecological states",
                  "EMP/SAP metabolic potential",
                  "TEA/OAP indices and EMP comparison",
                  "Diagnostic figures",
                  "Bryant figure set",
                  "kME, state scores, random forest",
                  "FuzzyForest feature selection",
                  "Global network topology",
                  "Integrated driver tiers",
                  "Climate sensitivity for candidate drivers",
                  "State-active subnetworks and bridge taxa",
                  "Functional driver master table",
                  "State transition meta-network figure",
                  "Driver quadrant figure",
                  "State bridge functional figure",
                  "Final narrative timeline")
)

write_header <- function(title) {
  cat("## ", title, "\n\n", sep = "")
}

cat_to_file <- function(expr) {
  sink(output_file)
  on.exit(sink(), add = TRUE)
  force(expr)
}

cat_to_file({
  cat("# ROCS Data and Output Summary\n\n")
  cat("- Generated by: `scripts/inspect_data.R`\n")
  cat("- Generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  cat("- Working directory: `", here(), "`\n", sep = "")
  cat("- Purpose: fast review of current output shapes, key counts, and known interpretation flags.\n\n")

  write_header("Current Pipeline Order")
  cat("`run_all.sh` should run these numbered steps. Numeric starts are normalized, so `--start 7`, `--start 07`, `--start 7b`, `--start 07b`, and `--start 010` are accepted where the normalized step exists.\n\n")
  md_table(pipeline_steps)

  write_header("Input Availability")
  input_status <- data.table(
    input = c("UPSTREAM$tax_damage", "UPSTREAM$metadata", "UPSTREAM$kegg_mods",
              "CLASS$prokaryote_function", "CLASS$thermo",
              "results/stage1/prokaryotes_vst.rds"),
    path = c(UPSTREAM$tax_damage, UPSTREAM$metadata, UPSTREAM$kegg_mods,
             CLASS$prokaryote_function, CLASS$thermo,
             file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
  )
  input_status[, status := vapply(path, exists_mark, character(1))]
  md_table(input_status)

  write_header("Stage 1 Sample Metadata")
  stage_meta_path <- file.path(RESULTS$stage1, "sample_metadata_stage1.tsv")
  stage_meta <- safe_fread(stage_meta_path)
  if (!is.null(stage_meta)) {
    cat("- Rows:", nrow(stage_meta), "\n")
    cat("- Columns:", ncol(stage_meta), "\n")
    cat("- Columns present:", paste(names(stage_meta), collapse = ", "), "\n")
    if (all(c("core", "y_bp") %in% names(stage_meta))) {
      core_summary <- stage_meta[, .(
        samples = .N,
        min_age_kyr = round(min(y_bp, na.rm = TRUE) / 1000, 3),
        max_age_kyr = round(max(y_bp, na.rm = TRUE) / 1000, 3)
      ), by = core][order(core)]
      cat("\n")
      md_table(core_summary)
    }
  } else {
    cat("- Missing:", stage_meta_path, "\n\n")
  }

  write_header("Prokaryote CLR/VST Matrix")
  vst_path <- file.path(RESULTS$stage1, "prokaryotes_vst.rds")
  vst <- if (file.exists(vst_path)) readRDS(vst_path) else NULL
  if (!is.null(vst)) {
    cat("- Dimensions:", paste(dim(vst), collapse = " x "), "(samples x taxa)\n")
    cat("- Range:", fmt_num(min(vst, na.rm = TRUE)), "to", fmt_num(max(vst, na.rm = TRUE)), "\n")
    cat("- Any NAs:", any(is.na(vst)), "\n")
    cat("- Zero-variance taxa:", sum(apply(vst, 2, var, na.rm = TRUE) == 0), "\n\n")
  } else {
    cat("- Missing:", vst_path, "\n\n")
  }

  write_header("WGCNA Outputs")
  mods <- safe_fread(file.path(RESULTS$stage1, "wgcna", "module_assignments.tsv"))
  if (!is.null(mods)) {
    cat("- Module-assignment rows:", nrow(mods), "\n")
    cat("- Non-grey modules:", mods[module != "grey", uniqueN(module)], "\n")
    cat("- Module distribution:\n\n")
    md_table(mods[, .N, by = module][order(-N)])
    if (!is.null(vst)) {
      cat("- VST taxa present in modules:", length(intersect(colnames(vst), mods$taxon)), "\n")
      cat("- VST taxa missing module assignment:", length(setdiff(colnames(vst), mods$taxon)), "\n\n")
    }
  }
  soft_power <- safe_fread(file.path(RESULTS$stage1, "wgcna", "soft_power.tsv"))
  if (!is.null(soft_power)) {
    cat("- Soft-power table rows:", nrow(soft_power), "\n")
    cat("- Soft-power columns:", paste(names(soft_power), collapse = ", "), "\n\n")
  }
  preservation <- safe_fread(file.path(RESULTS$stage1, "wgcna", "preservation.tsv"))
  if (!is.null(preservation)) {
    cat("- Preservation summary:\n\n")
    if ("preserved" %in% names(preservation)) {
      md_table(preservation[, .N, by = preserved][order(preserved)])
    }
    md_table(head(preservation[order(-Zsummary)], 10))
  }

  write_header("HMM State Outputs")
  states <- safe_fread(file.path(RESULTS$hmm, "hmm_states.tsv"))
  if (!is.null(states)) {
    cat("- State rows:", nrow(states), "\n")
    cat("- Numeric states:", states[, uniqueN(state)], "\n")
    cat("- Labels:", states[, uniqueN(label)], "\n\n")
    cat("Numeric state distribution:\n\n")
    md_table(states[, .N, by = state][order(state)])
    cat("Label distribution:\n\n")
    md_table(states[, .N, by = label][order(label)])
    dup_labels <- states[, .(numeric_states = paste(sort(unique(state)), collapse = ", "),
                             n_numeric_states = uniqueN(state)), by = label][n_numeric_states > 1]
    if (nrow(dup_labels) > 0) {
      cat("Label collision warning: these labels map to multiple numeric states.\n\n")
      md_table(dup_labels)
    }
  }
  bic <- safe_fread(file.path(RESULTS$hmm, "bic_comparison.tsv"))
  if (!is.null(bic)) {
    cat("BIC comparison:\n\n")
    md_table(bic)
  }
  hmm_valid <- safe_fread(file.path(RESULTS$hmm, "hmm_validation_metrics.tsv"))
  if (!is.null(hmm_valid)) {
    cat("Held-out validation metrics (GeoB25202_R2):\n\n")
    md_table(hmm_valid)
  }
  hmm_select <- safe_fread(file.path(RESULTS$hmm, "hmm_model_selection.tsv"))
  if (!is.null(hmm_select)) {
    cat("Hybrid model selection table:\n\n")
    md_table(hmm_select)
  }
  fingerprints <- safe_fread(file.path(RESULTS$hmm, "state_fingerprints.tsv"))
  if (!is.null(fingerprints)) {
    cat("State fingerprints:\n\n")
    md_table(fingerprints)
  }
  mis_cross <- safe_fread(file.path(RESULTS$hmm, "state_mis_crosstab.tsv"))
  if (!is.null(mis_cross)) {
    cat("State x glacial/interglacial crosstab:\n\n")
    md_table(mis_cross)
  }

  write_header("EMP/SAP and TEA Outputs")
  emp_sample <- safe_fread(file.path(RESULTS$emp, "emp_sap_per_sample.tsv"))
  emp_summary <- safe_fread(file.path(RESULTS$emp, "emp_sap_summary.tsv"))
  taxon_capacity <- safe_fread(file.path(RESULTS$emp, "taxon_dg_capacity.tsv"))
  if (!is.null(emp_sample)) {
    cat("- EMP/SAP sample rows:", nrow(emp_sample), "\n")
    cat("- EMP/SAP columns:", paste(names(emp_sample), collapse = ", "), "\n\n")
  }
  if (!is.null(emp_summary)) {
    cat("EMP/SAP summary:\n\n")
    md_table(emp_summary)
  }
  if (!is.null(taxon_capacity)) {
    cat("Top taxon thermodynamic capacities:\n\n")
    md_table(head(taxon_capacity[order(-total_dg_capacity)], 10))
  }
  tea_sample <- safe_fread(file.path(RESULTS$tea, "tea_indices_per_sample.tsv"))
  tea_climate <- safe_fread(file.path(RESULTS$tea, "tea_climate_models.tsv"))
  oap_taxon <- safe_fread(file.path(RESULTS$tea, "oap_per_taxon.tsv"))
  if (!is.null(tea_sample)) {
    cat("- TEA sample rows:", nrow(tea_sample), "\n")
    cat("- TEA columns:", paste(names(tea_sample), collapse = ", "), "\n\n")
  }
  if (!is.null(tea_climate)) {
    cat("TEA/EMP/SAP climate models:\n\n")
    md_table(tea_climate)
  }
  if (!is.null(oap_taxon) && "dominant_class" %in% names(oap_taxon)) {
    cat("OAP dominant class counts:\n\n")
    md_table(oap_taxon[, .N, by = dominant_class][order(-N)])
  }

  write_header("Taxon Importance Outputs")
  rf_perf <- safe_fread(file.path(RESULTS$importance, "rf_model_performance.tsv"))
  rf_imp <- safe_fread(file.path(RESULTS$importance, "rf_variable_importance.tsv"))
  ff_perf <- safe_fread(file.path(RESULTS$importance_fuzzy, "ff_model_performance.tsv"))
  ff_imp <- safe_fread(file.path(RESULTS$importance_fuzzy, "ff_variable_importance.tsv"))
  state_scores <- safe_fread(file.path(RESULTS$importance, "state_importance_scores.tsv"))
  kme <- safe_fread(file.path(RESULTS$importance, "taxon_kme.tsv"))
  if (!is.null(kme)) cat("- kME rows/cols:", nrow(kme), "x", ncol(kme), "\n")
  if (!is.null(state_scores)) cat("- State score rows/cols:", nrow(state_scores), "x", ncol(state_scores), "\n\n")
  if (!is.null(rf_perf)) {
    cat("Random forest performance:\n\n")
    md_table(rf_perf)
  }
  if (!is.null(rf_imp)) {
    cat("Top random forest taxa:\n\n")
    top_rf_col <- intersect(c("Overall", "importance", "variable_importance"), names(rf_imp))[1]
    if (!is.na(top_rf_col)) setorderv(rf_imp, top_rf_col, order = -1)
    md_table(head(rf_imp[, intersect(c("taxon", "Overall", "module"), names(rf_imp)), with = FALSE], 10))
  }
  if (!is.null(ff_perf)) {
    cat("FuzzyForest performance:\n\n")
    md_table(ff_perf)
  }
  if (!is.null(ff_imp)) {
    cat("Top FuzzyForest taxa:\n\n")
    md_table(head(ff_imp[order(-variable_importance),
                         intersect(c("taxon", "variable_importance", "module"), names(ff_imp)),
                         with = FALSE], 10))
  }

  write_header("Network and State-Network Outputs")
  net_stats <- safe_fread(file.path(RESULTS$network_stats, "network_metrics_summary.tsv"))
  if (!is.null(net_stats)) {
    cat("- Network metric rows:", nrow(net_stats), "\n")
    cat("- Network metric columns:", paste(names(net_stats), collapse = ", "), "\n\n")
    cat("Importance-type distribution:\n\n")
    md_table(net_stats[, .N, by = importance_type][order(-N)])
    cat("Z-P role distribution:\n\n")
    md_table(net_stats[, .N, by = role][order(-N)])
    cat("Top keystone candidates:\n\n")
    md_table(head(net_stats[is_potential_keystone == TRUE][order(rank_sum),
                            intersect(c("taxon", "module", "functional_group", "species",
                                        "rank_sum", "pagerank", "betweenness", "vulnerability"),
                                      names(net_stats)), with = FALSE], 10))
    cat("Hidden gems:\n\n")
    md_table(net_stats[is_hidden_gem == TRUE,
                       intersect(c("taxon", "module", "functional_group", "species", "pagerank",
                                   "z_score", "p_score"), names(net_stats)), with = FALSE])
  }
  state_net <- safe_fread(file.path(RESULTS$network_stats, "state_network_stats.tsv"))
  if (!is.null(state_net)) {
    cat("State-network metrics:\n\n")
    md_table(state_net)
  }
  bridge_taxa <- safe_fread(file.path(RESULTS$network_stats, "bridge_taxa_by_state.tsv"))
  if (!is.null(bridge_taxa)) {
    cat("Bridge taxa rows:", nrow(bridge_taxa), "\n")
    cat("Top bridge taxa per state:\n\n")
    top_bridges <- bridge_taxa[order(state, -inter_weight), head(.SD, 5), by = state]
    md_table(top_bridges[, intersect(c("state", "taxon", "species", "module", "driver_tier",
                                       "inter_edges", "inter_weight"), names(top_bridges)), with = FALSE])
  }

  write_header("Integrated Driver and Functional Linkage Outputs")
  drivers <- safe_fread(file.path(RESULTS$importance, "integrated_driver_summary.tsv"))
  master <- safe_fread(file.path(RESULTS$importance, "functional_driver_master.tsv"))
  climate <- safe_fread(file.path(RESULTS$importance, "climate_sensitivity_results.tsv"))
  state_func <- safe_fread(file.path(RESULTS$importance, "state_functional_enrichment.tsv"))
  if (!is.null(drivers)) {
    cat("- Integrated driver rows:", nrow(drivers), "\n")
    cat("- Integrated driver columns:", paste(names(drivers), collapse = ", "), "\n\n")
    cat("Driver-tier distribution:\n\n")
    md_table(drivers[, .N, by = driver_tier][order(driver_tier)])
    cat("Super-drivers:\n\n")
    md_table(drivers[is_super_driver == TRUE][order(-integrated_score),
                     intersect(c("taxon", "module", "functional_group", "species",
                                 "variable_importance", "topo_composite_percentile",
                                 "integrated_score"), names(drivers)), with = FALSE])
  }
  if (!is.null(climate)) {
    sig_n <- if ("fdr_d18O" %in% names(climate)) climate[!is.na(fdr_d18O) & fdr_d18O < PARAMS$fdr_threshold, .N] else NA_integer_
    cat("- Climate sensitivity rows:", nrow(climate), "\n")
    cat("- Climate-sensitive taxa at FDR <", PARAMS$fdr_threshold, ":", sig_n, "\n\n")
    md_table(head(climate[order(fdr_d18O)], 10))
  }
  if (!is.null(master)) {
    cat("- Functional driver master rows:", nrow(master), "\n")
    if ("is_functional_hub" %in% names(master)) {
      cat("- Functional hubs:", master[is_functional_hub == TRUE, .N], "\n\n")
    }
  }
  if (!is.null(state_func)) {
    cat("State functional enrichment:\n\n")
    md_table(state_func)
  }

  write_header("Output Group Inventory")
  inv <- results_inventory()
  if (nrow(inv) > 0) {
    cat("- Total files under `results/`:", nrow(inv), "\n")
    cat("- Total size under `results/`:", fmt_bytes(sum(inv$bytes_num, na.rm = TRUE)), "\n")
    cat("- Large SourceTracker feature-table directories are summarized by group; individual feature tables are intentionally not listed here.\n\n")
    md_table(output_group_summary(inv), max_rows = 40)
  } else {
    cat("- No files found under `results/`.\n\n")
  }

  write_header("Primary Output Table Schemas")
  schema_paths <- c(
    "results/stage1/sample_metadata_stage1.tsv",
    "results/stage1/prokaryotes_taxa_metadata.tsv",
    "results/stage1/wgcna/module_assignments.tsv",
    "results/stage1/wgcna/module_eigengenes.tsv",
    "results/stage1/wgcna/soft_power.tsv",
    "results/stage1/wgcna/preservation.tsv",
    "results/hmm/hmm_states.tsv",
    "results/hmm/state_fingerprints.tsv",
    "results/hmm/state_mis_crosstab.tsv",
    "results/hmm/bic_comparison.tsv",
    "results/hmm/hmm_validation_metrics.tsv",
    "results/hmm/hmm_model_selection.tsv",
    "results/hmm/pca_loadings.tsv",
    "results/emp/emp_sap_per_sample.tsv",
    "results/emp/emp_sap_summary.tsv",
    "results/emp/taxon_dg_capacity.tsv",
    "results/emp/taxon_sap.tsv",
    "results/tea/tea_indices_per_sample.tsv",
    "results/tea/oap_per_taxon.tsv",
    "results/tea/tea_climate_models.tsv",
    "results/tea/tea_vs_emp_correlations.tsv",
    "results/importance/rf_model_performance.tsv",
    "results/importance/rf_variable_importance.tsv",
    "results/importance/taxon_kme.tsv",
    "results/importance/state_importance_scores.tsv",
    "results/importance_fuzzy/ff_variable_importance.tsv",
    "results/importance_fuzzy/ff_model_performance.tsv",
    "results/network_stats/network_metrics_summary.tsv",
    "results/network_stats/state_network_stats.tsv",
    "results/network_stats/bridge_taxa_by_state.tsv",
    "results/importance/integrated_driver_summary.tsv",
    "results/importance/climate_sensitivity_results.tsv",
    "results/importance/functional_driver_master.tsv",
    "results/importance/state_functional_enrichment.tsv"
  )
  schema <- rbindlist(lapply(schema_paths, function(p) schema_row(here(p), p)), fill = TRUE)
  md_table(schema)

  write_header("Figure Inventory")
  if (exists("inv") && nrow(inv) > 0) {
    fig_dt <- inv[ext %in% c("png", "pdf"),
                  .(file, ext, bytes = fmt_bytes(bytes_num))][order(file)]
    cat("- Figure files across `results/`:", nrow(fig_dt), "\n\n")
    md_table(fig_dt, max_rows = 120)
  } else {
    cat("- No result inventory available for figure listing.\n\n")
  }

  write_header("Module Count and Color Divergence")
  current_modules <- if (!is.null(mods)) sort(unique(mods$module)) else character()
  current_non_grey <- setdiff(current_modules, "grey")
  cat("- Current output module colors:", paste(current_modules, collapse = ", "), "\n")
  cat("- Current non-grey module count:", length(current_non_grey), "\n")
  cat("- Current non-grey module colors:", paste(current_non_grey, collapse = ", "), "\n")
  cat("- Current output contains red module:", "red" %in% current_modules, "\n\n")

  module_classification <- data.table(
    file = c(
      "results/stage1/wgcna/module_assignments.tsv",
      "DATA_SUMMARY.md",
      "CODEX.SUMMARY.md",
      "README.md",
      "BryantsNotes.md",
      "CODEX.md",
      "GEMINI.md",
      "WORKFLOW_DOCUMENTATION.md",
      "old/13_fig_module_composition.R",
      "scripts/06b_bryantfigures.R",
      "scripts/13_state_transition_network.R",
      "scripts/16_final_story_visualization.R",
      "results/importance/taxon_kme.tsv"
    ),
    classification = c(
      "current 5 non-grey modules",
      "current 5-module memory file",
      "current 5-module memory file",
      "no explicit module count",
      "current hybrid-HMM notes file",
      "current 5-module + hybrid-HMM prose",
      "stale 6-module prose",
      "current 5-module + hybrid-HMM prose",
      "older 6 climate-linked / 7 total module script",
      "stale red-module palette only",
      "stale six-module/red layout",
      "stale red-module palette only",
      "stale downstream shape"
    ),
    evidence = c(
      paste(current_non_grey, collapse = ", "),
      "Generated/updated summary should be regenerated from current outputs.",
      "Lists current module counts and flags 6-module prose as drift.",
      "Mentions WGCNA modules without asserting a count.",
      "Includes 2026-05-05 hybrid HMM decision and holdout metrics.",
      "States 5 stable non-grey modules and hybrid HMM selection workflow.",
      "Says 6 stable modules: Turquoise, Blue, Brown, Yellow, Green, Red.",
      "Workflow graph/table states 5 modules and hybrid HMM metrics.",
      "Defines climate-linked modules as blue, brown, turquoise, yellow, black, red.",
      "Module palettes include red although current output has no red module.",
      "Graph layout is turquoise-blue-brown-yellow-green-red-turquoise.",
      "Module palette includes red although current output has no red module.",
      "Header includes kMEred although current output has no red module."
    )
  )
  md_table(module_classification)

  divergence_files <- c(
    "CODEX.md",
    "GEMINI.md",
    "WORKFLOW_DOCUMENTATION.md",
    "README.md",
    "BryantsNotes.md",
    "CODEX.SUMMARY.md",
    "old/13_fig_module_composition.R",
    "scripts/06b_bryantfigures.R",
    "scripts/13_state_transition_network.R",
    "scripts/16_final_story_visualization.R"
  )
  divergence_pattern <- paste(
    c("6 stable", "6 Stable", "6 Consensus", "7 modules",
      "climate-linked modules", "Climate-linked modules",
      "Current WGCNA output has", "5 non-grey", "red module",
      "turquoise.*blue.*brown.*yellow.*green.*red",
      "turquoise.*blue.*brown.*yellow.*black.*red"),
    collapse = "|"
  )
  divergence_hits <- rbindlist(lapply(divergence_files, function(path) {
    full_path <- here(path)
    if (!file.exists(full_path)) {
      return(data.table(file = path, line = NA_integer_, excerpt = "missing"))
    }
    txt <- readLines(full_path, warn = FALSE)
    idx <- grep(divergence_pattern, txt, ignore.case = TRUE)
    if (!length(idx)) {
      return(data.table(file = path, line = NA_integer_, excerpt = "no explicit module-count/color claim found"))
    }
    data.table(file = path, line = idx, excerpt = trimws(txt[idx]))
  }), fill = TRUE)
  cat("Detailed grep hits for auditability:\n\n")
  md_table(divergence_hits, max_rows = 80)

  write_header("Known Current Interpretation Flags")
  cat("- Current WGCNA output has ", if (!is.null(mods)) mods[module != "grey", uniqueN(module)] else "NA",
      " non-grey modules. Current module colors are ",
      if (!is.null(mods)) paste(sort(unique(mods$module)), collapse = ", ") else "NA", ".\n", sep = "")
  if (!is.null(states)) {
    cat("- Current HMM output has ", states[, uniqueN(state)], " numeric states and ",
        states[, uniqueN(label)], " labels.\n", sep = "")
  }
  if (!is.null(hmm_select) && nrow(hmm_select) > 0) {
    best_row <- hmm_select[order(-bic_ambiguous, -validation_logLik_per_sample, BIC)][1]
    cat("- Hybrid HMM metrics are available in `results/hmm/hmm_validation_metrics.tsv` and `results/hmm/hmm_model_selection.tsv`; current selected K is ",
        best_row$K, ".\n", sep = "")
  }
  cat("- `11_state_networks.R` filters active taxa from a global adjacency; treat these as state-active subnetworks, not independently inferred state-specific networks.\n")
  cat("- `run_all.sh` now includes the numbered driver/story scripts through step 16 and normalizes numeric starts plus suffix starts such as `7b`/`07b`; execution still requires a valid `Rscript` in the runtime environment.\n")
})

summary_lines <- readLines(output_file, warn = FALSE)
writeLines(sub("[[:space:]]+$", "", summary_lines), output_file, useBytes = TRUE)

message("Summary written to ", output_file)
