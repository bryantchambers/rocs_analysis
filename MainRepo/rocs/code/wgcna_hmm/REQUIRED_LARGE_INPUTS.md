# Required Large Inputs Not Tracked In GitHub

The `wgcna_hmm` workflow expects these files from `00_config.R`:

- `results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz`
- `data/functional/kegg-modules-summary-rocs.tsv.gz`

These two files are each larger than GitHub's 100 MB file limit, so they are not committed to this repository.

To run the full workflow, copy both files into the exact paths above after cloning.

You can also set explicit paths with environment variables:

- `WGCNA_HMM_INPUT_TAX_DAMAGE`
- `WGCNA_HMM_INPUT_KEGG_MODS`
- `WGCNA_HMM_INPUT_METADATA`
- `WGCNA_HMM_INPUT_XRF`
- `WGCNA_HMM_INPUT_PROK_FUNCTION`

`00_config.R` also includes a fallback to `/projects/caeg/people/ngm902/apps/repos/rocs/...` so collaborators on the same server can use your shared files without copying them.
