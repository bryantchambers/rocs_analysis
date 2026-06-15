# CLR Update Merge Audit

- Merge branch: `CLR-update`
- Source workflow root: `/src`
- Target workflow root: `/src/MainRepo/rocs`
- Canonical M git repo: `/src/MainRepo/rocs/.git`
- Legacy git metadata present but not used for this merge: `/src/MainRepo/rocs/.MAT-GIT`

## Git Note

`/src/MainRepo/rocs` currently contains both `.git` and `.MAT-GIT`.

- `.git` is the active operational repo for this merge and for future merge requests.
- `.MAT-GIT` is legacy metadata from an older merge workflow and is retained only as historical context.

All branch creation, validation, and merge-request work for this update should use `.git`.

## Merge Boundary

This merge updates M's workflow only through the WGCNA handoff.

- Bryant workflow remains untouched in `/src`.
- M's HMM and downstream workflow remain the default after WGCNA.
- Balanced sampling still applies to WGCNA module construction only; projected module eigengenes continue to feed M's HMM stage.

## Production WGCNA Input Contract

The active production input contract is:

- abundance feedstock: aggregated `tax_abund_tad`
- prevalence basis: `tax_abund_tad > 0`
- centering: sample-wise CLR with pseudocount `0.5`

Legacy filenames such as `clr_matrix_train_centered.rds` are preserved only for downstream path compatibility.

## Default and Preserved Profiles

The default balanced profile for this CLR update is:

|profile|power|deepSplit|mergeCutHeight|minModuleSize|observed neighborhood|
|---|---:|---:|---:|---:|---|
|`balanced_clr_default`|3|4|0.02|8|7 non-grey modules, `grey_pct ~ 48.55`|

Preserved balanced alternates:

|profile|power|deepSplit|mergeCutHeight|minModuleSize|observed neighborhood|
|---|---:|---:|---:|---:|---|
|`balanced_clr_alt001`|3|4|0.01|8|tied 7-module neighborhood|
|`balanced_clr_alt005`|3|4|0.05|8|tied 7-module neighborhood|
|`balanced_clr_fallback_p2_ds3`|2|3|0.02|6|6 non-grey modules, `grey_pct ~ 55.49`|
|`balanced_clr_fallback_p2_ds4`|2|4|0.02|6|6 non-grey modules, `grey_pct ~ 55.49`|

Historical profiles preserved for compatibility:

|profile|power|deepSplit|mergeCutHeight|minModuleSize|role|
|---|---:|---:|---:|---:|---|
|`original_clr_default`|2|4|0.02|6|runnable original-selection fallback under corrected CLR|
|`balanced_top3`|12|1|0.25|30|older balanced merge default|
|`original_exp3`|12|3|0.25|20|historical original-method profile; preserved for comparison|

Use the runnable original-selection fallback with:

```bash
WGCNA_HMM_INPUT_STRATEGY=original WGCNA_HMM_WGCNA_PROFILE=original_clr_default
```

Use the historical original exp3 profile with:

```bash
WGCNA_HMM_INPUT_STRATEGY=original WGCNA_HMM_WGCNA_PROFILE=original_exp3
```

## Workflow Contract

The merged workflow keeps:

- `code/wgcna_hmm/01_data_prep.R`: builds the corrected CLR matrix and WGCNA training manifest
- `code/wgcna_hmm/02_wgcna_main.R`: fits WGCNA on the selected balanced manifest by default
- `code/wgcna_hmm/02b_wgcna_stability.R`: runs bootstrap module stability diagnostics
- `code/wgcna_hmm/run_workflow.R`: runs WGCNA, stability, and M's downstream HMM/TEA/report stages

Supporting QC directories remain in repo:

- `InputQC/`
- `networkQC/`

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

Balanced-specific outputs remain:

- `results/wgcna_hmm/main/wgcna_training_samples.tsv`
- `results/wgcna_hmm/main/balanced_baseline_samples.tsv`
- `results/wgcna_hmm/main/balance_bin_quotas.tsv`
- `results/wgcna_hmm/main/balance_bin_availability.tsv`
- `results/wgcna_hmm/main/balance_excluded_samples.tsv`
- `results/wgcna_hmm/main/balance_design_summary.tsv`
- `results/wgcna_hmm/main/wgcna_stability/module_stability_bootstrap.tsv`
- `results/wgcna_hmm/main/wgcna_stability/module_stability_summary.tsv`
- `results/wgcna_hmm/main/wgcna_stability/module_size_sensitivity.tsv`
- `results/wgcna_hmm/main/wgcna_stability/WGCNA_STABILITY_REPORT.md`

## Run Commands

Development/build run with the CLR default:

```bash
WGCNA_HMM_INPUT_STRATEGY=balanced \
WGCNA_HMM_WGCNA_PROFILE=balanced_clr_default \
Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

Final run with the CLR default:

```bash
WGCNA_HMM_INPUT_STRATEGY=balanced \
WGCNA_HMM_WGCNA_PROFILE=balanced_clr_default \
Rscript code/wgcna_hmm/run_workflow.R --mode=final
```

Original fallback smoke test:

```bash
WGCNA_HMM_INPUT_STRATEGY=original \
WGCNA_HMM_WGCNA_PROFILE=original_clr_default \
Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

Balanced alternate smoke test:

```bash
WGCNA_HMM_INPUT_STRATEGY=balanced \
WGCNA_HMM_WGCNA_PROFILE=balanced_clr_alt001 \
Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

## Validation Checklist

- New default profile completes end-to-end in build mode.
- New default profile completes end-to-end in final mode.
- Original fallback still completes in build mode.
- Reports and run metadata record `balanced_clr_default` plus the corrected input contract.
- M downstream output filenames remain unchanged.

## Recovery Note

Recovery for this merge should use `/src/MainRepo/rocs/.git`, not `.MAT-GIT`.
