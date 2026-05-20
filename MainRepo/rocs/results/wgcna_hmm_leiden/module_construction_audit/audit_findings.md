# Hybrid WGCNA-consensus + Leiden module-construction audit

Generated: 2026-05-06

## Scope and assumptions

- Scope: audit and improve the Leiden-based hybrid module construction under `results/wgcna_hmm_leiden/`.
- Baseline `results/wgcna_hmm/` was left untouched.
- Assumption: `WGCNA::consensusDissTOMandTree()` returns a **dissimilarity-like** matrix in `consensusTOM` for this usage context (empirically supported here by distribution diagnostics).

## Key pathology in the original run

Original pattern (resolution 1.6):

- 25 non-grey modules
- 24 modules had exactly 24 taxa
- grey = 1076 taxa

This is pathological and strongly suggests graph-construction artifacts rather than biological structure.

## Audit findings

### 1) Mutual kNN vs union kNN: implementation/provenance mismatch

- Code used:
  - directed kNN retention, then `pmax(sparse_adj, t(sparse_adj))`
  - this is **union kNN**, not mutual kNN.
- Metadata/log text claimed “mutual kNN”.

Conclusion: **No algorithmic bug in symmetrization operation itself**, but a **provenance/labeling bug** (misreported graph policy).

### 2) `consensusDissTOMandTree` output interpretation was wrong for Leiden graph weights

- The hybrid code directly used `cons$consensusTOM` as similarity edge weights.
- Empirical check of off-diagonal values showed matrix concentrated near 1 (median ~0.99999), which is consistent with a dissimilarity-like matrix in this context.
- Using it as similarity causes nearest-neighbor selection to prefer **most dissimilar** pairs.

Conclusion: **Concrete algorithmic bug**: dissimilarity interpreted as similarity.

### 3) Why repeated 24-size modules happened

- With k=25 union kNN and min module size=20, at high resolution the partition produced many tiny communities.
- After filtering to keep only communities >=20, surviving communities clustered at ~24 taxa.
- This is a graph/design artifact amplified by the interpretation bug above (and by high resolution), not a meaningful module structure.

Conclusion: repeated 24s were due to **both**:

1. a real bug (dissimilarity treated as similarity), and
2. design pressure from sparse kNN + minimum-size filtering.

### 4) Resolution sweep/selection rule

- The old sweep did not itself create the bug, but it could select a pathological configuration because it operated on a malformed graph.
- Under corrected graph construction, sweep behavior changed substantially and selected a low-module solution (`resolution=1.0`, 6 modules, low grey).

Conclusion: selection logic was not the root bug, but it can mask upstream graph-construction problems.

### 5) Metadata/provenance mismatches identified

- “mutual kNN” was reported while union kNN was used.
- Graph representation text did not capture matrix interpretation mode.

Both are now fixed in metadata outputs.

## Conservative fixes implemented

1. Added explicit graph symmetrization parameter (`union` vs `mutual`) and corrected labeling.
2. Added explicit consensus matrix interpretation (`similarity`, `dissimilarity`, `auto`), with auto-heuristic and provenance capture.
3. Added optional eigengene-correlation post-merge hook (off by default), keeping change conservative.
4. Propagated new settings and provenance to `02_wgcna_main` and `10_tune_leiden_resolution` metadata.

## Focused exploration summary

See `module_construction_exploration.tsv` (24 settings over a small grid):

- `k_neighbors`: 25, 50, 100
- symmetrization: union, mutual
- min module size: 20, 30
- post-merge: none, corr_0.95
- fixed resolution: 1.6

Main observations:

- Corrected dissimilarity interpretation removed the exact-k-size pathology.
- Union kNN gave lower grey burden and more stable multi-module solutions than mutual kNN here.
- Larger k (50/100) tended to increase largest-module dominance; k=25 was most interpretable.
- Optional post-merge often made modules too coarse; not recommended as default.

## Recommended setting

See `recommended_setting.tsv`.

Recommended operational setting:

- consensus interpretation: `dissimilarity`
- `k_neighbors = 25`
- graph symmetrization: `union`
- `min_module_size = 20`
- post-merge: `none`
- `resolution = 1.0` (from corrected sweep)

## Rerun outcome at recommended setting

`run_workflow_leiden.R` rerun with corrected settings completed successfully.

Key resulting module sizes:

- grey: 148
- non-grey Leiden modules: 6 (`82, 141, 298, 342, 361, 425`)

This is substantially less pathological and more interpretable than the original 25 modules with repeated 24-size artifacts.

## Downstream compatibility checks

- Module names remain Leiden-derived (`leiden1..`, plus `grey`) in `main/module_assignments.tsv`.
- Eigengene columns remain `MEleiden*` in `main/module_eigengenes_main.tsv`.
- Main downstream stages (HMM/XRF/TEA/report/compare) completed in the rerun.

## Remaining caveats

- Consensus matrix interpretation relies on function behavior + empirical distribution; if upstream WGCNA semantics change, this should be revalidated.
- Very low grey burden can indicate broad modules; interpretability should be judged alongside downstream signal stability, not module count alone.
- Optional eigengene-merge remains available but should be used cautiously to avoid over-merging.
