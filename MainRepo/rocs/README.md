# ROCS collaboration README

This repository was prepared to share the `wgcna_hmm` workflow and its outputs with collaborators.

## What is included in this repo

- Workflow scripts under `code/wgcna_hmm/`
- Supporting scripts under other `code/` subfolders
- `wgcna_hmm` outputs:
  - `results/wgcna_hmm/`
  - `results/wgcna_hmm_leiden/`
- Required small/medium input files used by `wgcna_hmm`:
  - `data/metadata_v5.tsv`
  - `data/combined_xrf_not_normalized.tsv`
  - `results/microbial/wgcna/classification/prokaryote_function_assigned.tsv`

## What is intentionally NOT included

### Two required large input files are not in GitHub
GitHub rejects files larger than 100 MB. These required inputs exceed that limit:

- `data/functional/kegg-modules-summary-rocs.tsv.gz` (~210 MB)
- `results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz` (~138 MB)

## How this was solved

`code/wgcna_hmm/00_config.R` uses:

1. Relative paths for files that are in this repo
2. Absolute paths for the two large files that are not in GitHub:
   - `/projects/caeg/people/ngm902/apps/repos/rocs/data/functional/kegg-modules-summary-rocs.tsv.gz`
   - `/projects/caeg/people/ngm902/apps/repos/rocs/results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz`

## Running the workflow

From repo root:

```bash
Rscript code/wgcna_hmm/run_workflow.R --mode=build
```

For final WGCNA preservation/stability settings:

```bash
Rscript code/wgcna_hmm/run_workflow.R --mode=final
```

Leiden alternative workflow:

```bash
Rscript code/wgcna_hmm/run_workflow_leiden.R
```

## Bryant WGCNA/QC merge

This working branch includes Bryant's WGCNA and QC updates through the WGCNA handoff:

- `code/wgcna_hmm/01_data_prep.R` writes a WGCNA training manifest and balanced core-age-bin design tables.
- `code/wgcna_hmm/02_wgcna_main.R` now defaults to the balanced CLR profile `balanced_clr_default`.
- The production WGCNA input contract is aggregated `tax_abund_tad` with sample-wise CLR centering and prevalence on `tax_abund_tad > 0`.
- Alternate balanced CLR profiles are preserved for the tied `7`-module neighborhood and the nearby `6`-module fallback neighborhood.
- The original sample-selection path is preserved and now defaults to a CLR-era runnable fallback profile; the historical `original_exp3` profile remains available for comparison.
- `code/wgcna_hmm/02b_wgcna_stability.R` adds bootstrap module stability diagnostics before HMM.
- Top-level `InputQC/` and `networkQC/` contain the supporting QC analyses and reports.
- M's HMM and downstream scripts remain the default workflow after WGCNA.

Important HMM note: balanced sampling currently controls WGCNA module construction only. Module eigengenes are projected back to all main-window samples before M's HMM stage. HMM input balancing is expected to need review with M before final biological claims.

See `MERGE_AUDIT.md` for the merge boundary, file changes, output contract, and recovery notes.

Git note: `/src/MainRepo/rocs/.git` is the canonical repo for current merge work. `.MAT-GIT` is legacy metadata from an older merge workflow and is not the active repo for this branch.

## Notes

- Run scripts from repo root so relative paths resolve correctly.
- If running outside this shared server path, edit `code/wgcna_hmm/00_config.R` and update those two absolute paths.
