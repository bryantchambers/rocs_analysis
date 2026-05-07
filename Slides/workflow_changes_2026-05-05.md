# ROCS Workflow Update
## Pipeline Changes and Outcomes (2026-05-05)

- Scope: updates to WGCNA stability, hybrid HMM state discovery, and downstream biological interpretation scripts.
- Focus scripts:
  - `02_wgcna.R` (edited)
  - `02b_wgcna_stability.R` (new)
  - `03_hmm_states.R` (edited)
  - `07b_taxon_importance_fuzzy.R`, `08_network_statistics.R`, `09_driver_integration.R`
  - `10_climate_sensitivity.R`, `11_state_networks.R`, `12_functional_linkage.R`

---

# Executive Summary

## What Changed

1. WGCNA now has mode-aware validation depth (`build` vs `final`) and better validation outputs.
2. Added dedicated WGCNA stability script with bootstrap reruns + novice-readable report.
3. HMM moved to hybrid train/holdout/refit selection for stability on damaged aDNA.
4. Full pipeline rerun completed in `final` mode with updated state/module outputs.

---

# Pipeline Order Change

## `run_all.sh`

- Added new step:
  - `02b | WGCNA stability diagnostics | 02b_wgcna_stability.R`
- Added mode passthrough:
  - `bash run_all.sh --mode build|final`
  - mode is passed to steps `02` and `02b`.

---

# Config Updates

## `config.R` Additions

- New output directory:
  - `RESULTS$wgcna_stability = results/wgcna_stability`
- New WGCNA runtime controls:
  - `wgcna_run_mode` (`build` or `final`)
  - `wgcna_preservation_permutations_build = 200`
  - `wgcna_preservation_permutations_final = 700`
  - `wgcna_stability_bootstrap_build = 30`
  - `wgcna_stability_bootstrap_final = 120`
  - `wgcna_stability_age_grid_points = 100`

---

# `02_wgcna.R` Overview

## Core Process

1. Load CLR/VST matrix + sample metadata.
2. Select soft power from signed-R2 criteria across training cores.
3. Build consensus modules on `ST8`, `ST13`, `GeoB25202_R1`.
4. Compute eigengenes for train + validation (`GeoB25202_R2`).
5. Run module preservation (`R1 -> R2`).
6. Save module and eigengene artifacts.

---

# `02_wgcna.R` New Settings and Logic

## Changes Implemented

- New CLI mode argument:
  - `Rscript scripts/02_wgcna.R --mode=build|final`
- Mode-aware permutation depth:
  - build: `nPermutations = 200`
  - final: `nPermutations = 700`
- Preservation output split:
  - technical modules (`grey`, `gold`)
  - biological modules (non-grey/gold)
- Added age-aligned eigengene concordance:
  - interpolation on common age grid
  - Pearson, Spearman, RMSE metrics.

---

# `02_wgcna.R` Output Structure

## Key Outputs

- Existing:
  - `results/stage1/wgcna/module_assignments.tsv`
  - `results/stage1/wgcna/module_eigengenes.tsv`
  - `results/stage1/wgcna/soft_power.tsv`
  - `results/stage1/wgcna/preservation.tsv`
  - `results/stage1/wgcna/consensus_wgcna.rds`
- New:
  - `results/stage1/wgcna/preservation_biological.tsv`
  - `results/stage1/wgcna/eigengene_concordance_age_aligned.tsv`

---

# `02b_wgcna_stability.R` Overview (New)

## Purpose

- Quantify module stability under resampling.
- Summarize preservation and cross-core eigengene concordance.
- Produce readable report for method QA before biological claims.

## Core Steps

1. Load baseline modules + WGCNA outputs.
2. Bootstrap reruns of consensus WGCNA (mode-aware count).
3. Match baseline modules to best-overlap bootstrap modules.
4. Compute Jaccard stability distributions.
5. Write tabular outputs + markdown report.

---

# `02b_wgcna_stability.R` Output Structure

## New Artifacts (`results/wgcna_stability/`)

- `module_stability_bootstrap.tsv`
- `module_stability_summary.tsv`
- `module_size_sensitivity.tsv`
- `WGCNA_STABILITY_REPORT.md`

## Read-first File

- `WGCNA_STABILITY_REPORT.md`:
  - interpretation guide,
  - module-level stability table,
  - preservation summary,
  - age-aligned concordance summary,
  - next-action guidance.

---

# `03_hmm_states.R` Overview

## Hybrid State Discovery (Edited)

1. Residualize eigengenes within core.
2. PCA fit on training cores only.
3. Fit candidate HMMs (`K=2..5`) on training cores.
4. Score held-out `GeoB25202_R2` with fixed parameters.
5. Select K from BIC-ambiguous models using holdout + stability.
6. Refit selected K on all cores for final state labels.

---

# `03_hmm_states.R` New Metrics and Outputs

## New Files

- `results/hmm/hmm_validation_metrics.tsv`
- `results/hmm/hmm_model_selection.tsv`

## Key Selection Metrics

- train BIC and `delta_bic`
- holdout `validation_logLik_per_sample`
- `validation_mean_max_posterior`
- `validation_self_transition_rate`
- `validation_switches_per_100`

---

# Key Outcome: WGCNA + HMM

## Final-mode run outcomes

- WGCNA:
  - 5 non-grey modules (`blue`, `brown`, `green`, `turquoise`, `yellow`)
  - biological preservation: 4 strong, 1 moderate.
- HMM:
  - train-BIC best: `K=5`
  - hybrid-selected: `K=4`
  - final outputs now use stable 4-state solution.

---

# `07b_taxon_importance_fuzzy.R`

## Process Overview

1. Align taxa, modules, and HMM labels.
2. Train/test split for state prediction.
3. FuzzyForest feature selection by module structure.
4. Export model performance + top variable importance.

## Outputs

- `results/importance_fuzzy/ff_model_performance.tsv`
- `results/importance_fuzzy/ff_variable_importance.tsv`
- diagnostic plot and model object.

---

# `08_network_statistics.R`

## Process Overview

1. Build adjacency/TOM-based network from WGCNA-aligned taxa.
2. Compute Z-P roles and centralities (PageRank, closeness, betweenness).
3. Compute bridging centrality + nodal vulnerability.
4. Annotate taxa and classify influence types.

## Outputs

- `results/network_stats/network_metrics_summary.tsv`
- Z-P plots and state/network summary figures.

---

# `09_driver_integration.R`

## Process Overview

1. Merge predictive importance + topology + function.
2. Compute integrated scores and percentile composites.
3. Assign tiers:
  - Tier 1 super-driver,
  - Tier 2 high potential,
  - Tier 3 predictive specialist.

## Outputs

- `results/importance/integrated_driver_summary.tsv`
- super-driver plots and supporting summaries.

---

# `10_climate_sensitivity.R`

## Process Overview

1. Select candidate drivers from integrated summary.
2. Fit climate sensitivity models (GLS/CAR1 design).
3. Estimate `d18O` effects and FDR-adjusted significance.

## Outputs

- `results/importance/climate_sensitivity_results.tsv`
- volcano/sensitivity figures.

---

# `11_state_networks.R`

## Process Overview

1. Build state-active subnetworks from global network + state activity.
2. Quantify density/transitivity/inter-module ratio by state.
3. Identify bridge taxa within each state.

## Outputs

- `results/network_stats/state_network_stats.tsv`
- `results/network_stats/bridge_taxa_by_state.tsv`

---

# `12_functional_linkage.R`

## Process Overview

1. Join integrated drivers with EMP/TEA and climate sensitivity.
2. Identify functional super-hubs.
3. Summarize state bridge functional composition.

## Outputs

- `results/importance/functional_driver_master.tsv`
- `results/importance/state_functional_enrichment.tsv`
- functional linkage plots.

---

# Runtime and Operations

## Observed Full Final-mode Runtime

- Full pipeline (`01` through `16`, with `02b`) completed in ~53 minutes.
- Longest blocks:
  - `02` preservation permutations (700),
  - `02b` stability bootstraps (120).

## Recommended Operation

- Daytime iteration: `--mode build`
- End-of-day reproducibility/stability check: `--mode final`

---

# Risks and Follow-up

## Remaining Technical Notes

- Some figure scripts emit non-fatal Unicode/font warnings.
- LR04 input reader warns on one malformed line; run continues.
- Stability report should be reviewed before high-confidence biological interpretation.

## Next Step

- Add a one-page “decision gate” checklist:
  - minimum stability thresholds,
  - minimum concordance requirements,
  - holdout/HMM quality criteria.

