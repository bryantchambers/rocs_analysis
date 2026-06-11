# Bryant WGCNA/QC Merge Audit

- Merge branch: `bryant-balanced-wgcna-default`
- Source workflow root: `/src`
- Target workflow root: `/src/MainRepo/rocs`
- M git metadata: `.MAT-GIT`
- Baseline M commit before merge work: `891a8a4 Simplify wgcna_hmm input paths: relative for tracked files, absolute for large files`

## Merge Boundary

This merge intentionally updates M's workflow only through the WGCNA handoff.

- Bryant workflow remains untouched in `/src`.
- M's HMM and downstream workflow remain the default after WGCNA.
- Bryant's hybrid HMM K-selection/validation is not merged yet; it should be compared after this WGCNA/QC merge is stable.

## WGCNA Production Setting

The active WGCNA input contract is:

- abundance feedstock: aggregated `tax_abund_tad`
- prevalence basis: `tax_abund_tad > 0`
- centering: sample-wise CLR with pseudocount `0.5`

Legacy file names such as `clr_matrix_train_centered.rds` are preserved only for downstream path compatibility.

The default WGCNA input strategy is now `balanced`, using balancednetwork-selected `top3` parameters:

|parameter|value|
|---|---:|
|soft_power|12|
|deepSplit|1|
|mergeCutHeight|0.25|
|minModuleSize|30|

Rationale: balanced `top3` reduces ST8/core-age imbalance and had stronger balanced bootstrap stability while retaining reasonable preservation and kME.

The original networkQC-selected `exp3` setting is preserved as a fallback:

|parameter|value|
|---|---:|
|soft_power|12|
|deepSplit|3|
|mergeCutHeight|0.25|
|minModuleSize|20|

Use it with `WGCNA_HMM_INPUT_STRATEGY=original` or `WGCNA_HMM_WGCNA_PROFILE=original_exp3`.

## Critical HMM Note

Balanced sampling is currently applied to WGCNA module construction only.

M's HMM stage still consumes all projected main-window module eigengenes. This is intentional for this merge to preserve M's downstream workflow contract, but it is not a final methodological decision. HMM input balancing needs direct review with M before final biological claims about HMM state structure.

## Edited Files

- `.gitignore`
  - Tracks top-level `InputQC/` and `networkQC/` in M's repo.
- `code/wgcna_hmm/00_config.R`
  - Adds balanced WGCNA defaults and preserves original exp3 as a fallback profile.
  - Adds `WGCNA_HMM_INPUT_STRATEGY=balanced|original`.
  - Adds `WGCNA_HMM_WGCNA_PROFILE=balanced_top3|original_exp3|custom`.
  - Adds build/final mode controls.
  - Adds preservation/bootstrap/age-grid parameters.
  - Adds `DIRS$qc` and `DIRS$wgcna_stability`.
  - Honors documented input path environment variables.
- `code/wgcna_hmm/01_data_prep.R`
  - Builds WGCNA input from `tax_abund_tad` with sample-wise CLR.
  - Writes the WGCNA input metadata and training manifest.
  - Writes balanced core-age-bin QC tables when balanced mode is active.
- `code/wgcna_hmm/02_wgcna_main.R`
  - Fits WGCNA on the selected balanced training manifest by default.
  - Projects module eigengenes back to all main-window samples for M's HMM stage.
  - Uses build/final preservation permutation settings.
  - Adds biological-only preservation output.
  - Adds age-aligned R1/R2 eigengene concordance output.
  - Preserves M's training-basis validation eigengene projection and file names.
- `code/wgcna_hmm/02b_wgcna_stability.R`
  - New stage for bootstrap module stability diagnostics.
  - Resamples within balanced core-age-bin quotas in balanced mode.
  - Writes outputs to `results/wgcna_hmm/main/wgcna_stability/`.
- `code/wgcna_hmm/run_workflow.R`
  - Adds `--mode build|final`.
  - Runs `02b_wgcna_stability.R` after WGCNA and before HMM.
  - Writes `results/.../logs/workflow_progress.tsv` for run monitoring.

## Copied QC Directories

- `InputQC/`
  - Copied from Bryant workflow.
  - Contains input-depth, low-detection, rarefaction/depth, core-age imbalance, and ST8 low-taxa sensitivity analyses.
- `networkQC/`
  - Copied from Bryant workflow.
  - Contains WGCNA parameter sweep, full evaluation, kME/topology review, Leiden comparison, graph diagnostics, and reports.

## Main Output Contract Preserved

M downstream scripts should still consume:

- `results/wgcna_hmm/main/module_assignments.tsv`
- `results/wgcna_hmm/main/module_eigengenes_main.tsv`
- `results/wgcna_hmm/main/module_eigengenes_training.tsv`
- `results/wgcna_hmm/main/module_eigengenes_validation_projection.tsv`
- `results/wgcna_hmm/main/eigengene_projection_basis_main.tsv`
- `results/wgcna_hmm/main/module_preservation_validation.tsv`
- `results/wgcna_hmm/main/hmm_states_main.tsv`
- `results/wgcna_hmm/main/state_fingerprints_main.tsv`

New WGCNA/QC outputs include:

- `results/wgcna_hmm/main/wgcna_training_samples.tsv`
- `results/wgcna_hmm/main/balanced_baseline_samples.tsv`
- `results/wgcna_hmm/main/balance_bin_quotas.tsv`
- `results/wgcna_hmm/main/balance_bin_availability.tsv`
- `results/wgcna_hmm/main/balance_excluded_samples.tsv`
- `results/wgcna_hmm/main/balance_design_summary.tsv`
- `results/wgcna_hmm/main/module_preservation_validation_biological.tsv`
- `results/wgcna_hmm/main/eigengene_concordance_age_aligned.tsv`
- `results/wgcna_hmm/main/wgcna_stability/module_stability_bootstrap.tsv`
- `results/wgcna_hmm/main/wgcna_stability/module_stability_summary.tsv`
- `results/wgcna_hmm/main/wgcna_stability/module_size_sensitivity.tsv`
- `results/wgcna_hmm/main/wgcna_stability/WGCNA_STABILITY_REPORT.md`

## Run Commands

Development/build run:

```bash
WGCNA_HMM_INPUT_STRATEGY=balanced Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

Final run:

```bash
WGCNA_HMM_INPUT_STRATEGY=balanced Rscript code/wgcna_hmm/run_workflow.R --mode=final
```

Original exp3 fallback smoke test:

```bash
WGCNA_HMM_INPUT_STRATEGY=original WGCNA_HMM_WGCNA_PROFILE=original_exp3 Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

## Validation Runs

Balanced-default build validation:

```bash
WGCNA_HMM_RESULTS_SUFFIX=balanced_default_build_20260601_codex \
WGCNA_HMM_INPUT_TAX_DAMAGE=/src/results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz \
WGCNA_HMM_INPUT_KEGG_MODS=/src/data/functional/kegg-modules-summary-rocs.tsv.gz \
WGCNA_HMM_INPUT_STRATEGY=balanced \
WGCNA_HMM_WGCNA_PROFILE=balanced_top3 \
Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

Result: completed end-to-end through `06_report.R`.

Key checks:

- Main samples: `214`
- WGCNA selected training samples: `57`
- Selected training counts: `19` per training core
- Non-grey modules: `5`
- Module sizes: `turquoise=347`, `blue=312`, `brown=138`, `yellow=88`, `green=40`, `grey=872`
- HMM states: `5`
- HMM input balancing status: `pending_review_with_M`

Original fallback smoke validation:

```bash
WGCNA_HMM_RESULTS_SUFFIX=original_fallback_smoke_20260601_codex \
WGCNA_HMM_INPUT_TAX_DAMAGE=/src/results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz \
WGCNA_HMM_INPUT_KEGG_MODS=/src/data/functional/kegg-modules-summary-rocs.tsv.gz \
WGCNA_HMM_INPUT_STRATEGY=original \
WGCNA_HMM_WGCNA_PROFILE=original_exp3 \
WGCNA_HMM_PRESERVATION_PERMUTATIONS_BUILD=5 \
WGCNA_HMM_STABILITY_BOOTSTRAP_BUILD=2 \
Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

Result: completed end-to-end through `06_report.R`.

Key checks:

- Main samples: `214`
- Original WGCNA training samples: `189`
- Original training counts: `ST8=115`, `ST13=48`, `GeoB25202_R1=26`
- Non-grey modules: `8`
- HMM states: `5`

## Recovery Notes

- M's original `main` remains recoverable from `.MAT-GIT`.
- This merge branch can be inspected with:

```bash
git --git-dir=.MAT-GIT --work-tree=. status --short --branch
```

- Bryant-side source directories are not required at runtime for the merged M workflow; balanced design logic is integrated into `code/wgcna_hmm/`.
