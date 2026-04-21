#!/usr/bin/env Rscript
# config.R — centralized paths and parameters for rocs-repro
# Source this at the top of every script: source(here::here("config.R"))

library(here)

# ── Input data (upstream pipeline outputs, read-only) ────────────────────────

BASE <- "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs"

UPSTREAM <- list(
  # Per-taxon read counts, pre-filtered to damaged reads (is_dmg column)
  tax_damage = "/maps/projects/caeg/people/ngm902/apps/repos/rocs_marine_cores_old/results/microbial/taxonomy/dmg-summary-ssp_selected.tsv.gz",
  # Sample metadata (label, core, y_bp, mis, temp_complete)
  metadata   = file.path(BASE, "data/metadata_v4.txt"),
  # KEGG module completeness per genome (includes enzyme_hits_in_module)
  kegg_mods  = file.path(BASE, "external/functionalDB/kegg-modules-summary-27-03-2026.tsv.gz")
)

# ── Classification reference files ───────────────────────────────────────────

CLASS_DIR <- file.path(BASE, "analysis_wgcna/classification")

CLASS <- list(
  prokaryote_function = file.path(CLASS_DIR, "prokaryote_function_assigned.tsv"),
  eukaryote           = file.path(CLASS_DIR, "eukaryote_classification.tsv"),
  thermo              = file.path(CLASS_DIR, "dpi_reaction_thermodynamics.tsv"),
  kegg_module_class   = file.path(CLASS_DIR, "kegg_module_classification.tsv"),
  damage_tiers        = file.path(CLASS_DIR, "damage_authentication_tiers.tsv")
)

# ── Output directories ────────────────────────────────────────────────────────

# Pre-computed stage1 results from the original analysis (read-only)
OLD <- list(
  stage1   = "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1",
  wgcna    = "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1/wgcna_prokaryotes",
  hmm      = "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1/hmm_states",
  emp_sap  = "/maps/projects/caeg/people/kbd606/scratch/mateu-rocs/analysis_wgcna/results/stage1/emp_sap"
)



# New outputs (written here)
RESULTS <- list(
  stage1  = here("results", "stage1"),
  hmm     = here("results", "hmm"),
  emp     = here("results", "emp"),
  tea     = here("results", "tea"),
  figures = here("results", "figures")
)

for (d in RESULTS) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Analysis parameters ───────────────────────────────────────────────────────

PARAMS <- list(
  # Sample filtering
  stage1_max_age_kyr = 150,           # Stage-1 training window
  # All cores loaded in data prep; split in WGCNA
  all_cores          = c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"),
  stage1_cores       = c("ST8", "ST13", "GeoB25202_R1"),
  validation_core    = "GeoB25202_R2",
  excluded_samples   = "LV3003046968",

  # Damage filter (CCC criterion)
  dmg_rho_c  = 0.75,
  dmg_c_b    = 0.90,
  dmg_pval   = 0.10,
  dmg_euk_ratio = 0.50,

  # Data prep
  prevalence_min_samples = 10,    # absolute sample count (same as original)
  transform              = "clr", # CLR — original uses this, not VST

  # WGCNA
  wgcna_min_module_size  = 20,
  wgcna_deep_split       = 2,
  wgcna_merge_cut_height = 0.15,  # original value
  wgcna_soft_power       = NULL,  # auto-selected (scale-free R2 >= 0.80)
  wgcna_eigengene_r_min  = 0.70,
  wgcna_module_zscore    = 2,

  # HMM
  hmm_k       = 5,                # validated number of states
  hmm_n_iter  = 500,
  hmm_n_start = 10,

  # FlashWeave
  fw_sensitive    = FALSE,
  fw_heterogeneous = TRUE,
  fw_n_bootstrap  = 50,

  # EMP
  emp_default_dg = 100,           # default heterotrophy ΔG (kJ/mol)

  # TEA indices
  oap_completeness_threshold = 0.30,

  # Statistics
  fdr_threshold = 0.05,
  seed          = 42
)

# ── TEA module definitions (matching 14-tea-indices.R from sapropel project) ─

TEA <- list(
  # OAP modules: module → (TEA class, E0' mV, discrimination weight)
  oap_modules = list(
    M00155 = list(class = "O2",  eo_mv =  820, disc = 1.00),
    M00154 = list(class = "O2",  eo_mv =  820, disc = 1.00),
    M00416 = list(class = "O2",  eo_mv =  820, disc = 1.00),
    M00417 = list(class = "O2",  eo_mv =  820, disc = 1.00),
    M00156 = list(class = "O2",  eo_mv =  820, disc = 0.50),
    M00153 = list(class = "O2",  eo_mv =  820, disc = 0.25),
    M00529 = list(class = "NO3", eo_mv =  740, disc = 1.00),
    M00530 = list(class = "NO3", eo_mv =  360, disc = 1.00),
    M00973 = list(class = "ANX", eo_mv =  350, disc = 1.00),
    M00596 = list(class = "SO4", eo_mv = -217, disc = 1.00),
    M00567 = list(class = "CH4", eo_mv = -244, disc = 1.00),
    M00357 = list(class = "CH4", eo_mv = -244, disc = 1.00),
    M00356 = list(class = "CH4", eo_mv = -244, disc = 1.00),
    M00563 = list(class = "CH4", eo_mv = -244, disc = 1.00)
  ),
  # KOs for ratio-based indices
  target_kos = c(
    "K00399", "K00401", "K00402",  # mcrA/B/G — MGI
    "K14080",                      # pmoA     — MFI
    "K00370",                      # narG     — DCI
    "K02305",                      # norB     — DCI
    "K00376",                      # nosZ     — DCI
    "K00425",                      # cydA     — MII
    "K02274", "K02276",            # coxA/C   — MII
    "K00394", "K00395",            # aprA/B   — SRPI
    "K00958",                      # sat      — SRPI
    "K11180", "K11181"             # dsrA/B   — SRPI
  )
)

# ── Climate reference data ───────────────────────────────────────────────────

CLIMATE <- list(
  mis_boundaries = file.path(BASE, "results/network_global/mis_stage_boundaries.tsv"),
  lr04_stack     = "/projects/fernandezguerra/people/ngm902/ROCS/associated_data/Lisiecki2005_copy.txt"
)

# ── Figure palettes ──────────────────────────────────────────────────────────

PALETTES <- list(
  # MIS climate periods: glacial vs interglacial
  mis_climate = c(
    glacial      = "#6cd3f5",
    interglacial = "#fdb57ee8"
  ),
  # Prokaryote functional groups
  prok_function = c(
    Core_heterotrophy           = "#1B9E77",
    Methanogenesis_diagenetic   = "#D95F02",
    Pelagic_heterotroph_MGII    = "#7570B3",
    Nitrification_AOA           = "#E7298A",
    Unknown                     = "#CCCCCC"
  )
)