# Leiden resolution sweep selection summary

Generated: 2026-05-06 21:34:58

## Provenance
- Results root: `/maps/projects/caeg/people/ngm902/apps/repos/rocs/results/wgcna_hmm_leiden`
- Script: `code/wgcna_hmm/10_tune_leiden_resolution.R`
- Training cores: ST8, ST13, GeoB25202_R1
- Min module size threshold: 20
- Consensus graph: union kNN sparse graph from interpreted consensus matrix (k=25)
- Consensus matrix interpretation: dissimilarity
- Previous run metadata requested resolution: 1.6
- Previous run metadata used resolution in practice: 1.6

## Resolutions tested
- 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0

## Selection rule (conservative and explicit)
- 1) Candidate must pass guardrails: non-grey modules in [3,30], grey fraction <= 0.70, largest non-grey module <= 60% of taxa, HMM proxy runnable.
- 2) Among feasible candidates, retain those within 0.01 modularity of the best feasible modularity.
- 3) From that near-optimal set, choose the fewest non-grey modules (interpretability).
- 4) Tie-break by higher modularity, then lower resolution.

## Selected resolution
- Selected resolution: **1.000**
- Selection status: `selected_from_feasible`
- Modularity: 0.57014
- Non-grey modules: 6
- Grey fraction: 0.0824
- Non-grey size summary (min/median/max): 82 / 319.5 / 426

## Notes
- This is a sensitivity/alternative partition workflow; baseline results/wgcna_hmm remain unchanged.
- Full sweep table is in `leiden_resolution_sweep.tsv`.
