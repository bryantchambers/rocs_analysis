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
Rscript code/wgcna_hmm/run_workflow.R
```

Leiden alternative workflow:

```bash
Rscript code/wgcna_hmm/run_workflow_leiden.R
```

## Notes

- Run scripts from repo root so relative paths resolve correctly.
- If running outside this shared server path, edit `code/wgcna_hmm/00_config.R` and update those two absolute paths.
