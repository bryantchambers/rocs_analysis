# Original vs Balanced WGCNA Method Comparison

## Executive Summary

- Balanced best setting: `balanced_top3` (power=12, deepSplit=1, mergeCutHeight=0.25, minModuleSize=30).
- Original best setting: `exp3` (power=12, deepSplit=3, mergeCutHeight=0.25, minModuleSize=20).
- Balanced best has stronger bootstrap stability: 0.604 vs original exp3 0.400.
- Original exp3 has lower grey burden: 28.66% vs balanced best 48.53%.
- Original exp3 has stronger TOM separation: 3.532 vs balanced best 2.608.
- Current downstream HMM/driver/functional outputs represent the original/current pipeline, not a full balanced downstream rerun.

## Sampling Fairness

|method|training samples|ST8 fraction|core count CV|age-bin count CV mean|
|---|---:|---:|---:|---:|
|original|189|0.608|0.736|0.920|
|balanced|57|0.333|0.000|0.000|

## Key Method Comparison

|method|grey %|non-grey modules|bootstrap Jaccard|Pearson concordance|bio median kME|TOM separation|
|---|---:|---:|---:|---:|---:|---:|
|balanced_best_top3|48.53|5|0.604|0.835|0.706|2.608|
|original_best_exp3|28.66|8|0.400|0.857|0.703|3.532|
|original_baseline|66.56|5|0.351|0.847|0.783|48.011|

## Metric Winners

- Sampling fairness: balanced (`ST8 fraction 0.333`, core count CV 0.000).
- Bootstrap stability: balanced best (`0.604`).
- Grey burden and module recovery: original exp3 (`28.66%` grey, `8` non-grey modules).
- TOM separation among non-baseline candidates: original exp3 (`3.532`).
- Downstream biology availability: original/current pipeline (`9/9` artifacts available); balanced downstream rerun not yet available.

## Interpretation

- Balanced is more sampling-fair and, in the completed 500-bootstrap run, more label-stable than original exp3.
- Original exp3 is less grey, has more biological modules, stronger TOM separation, and currently owns the downstream biology evidence.
- The fair answer is not simply that one method wins. Balanced should be treated as the stronger sampling-control network; original exp3 remains the stronger current end-to-end biological network until balanced downstream HMM/driver analyses are rerun.
- If balanced downstream biology recovers the same state/driver story, promote balanced or use it as the primary sensitivity-supported network. If it loses major ecological signal, keep original exp3 with an explicit ST8-dominance caveat.

## Workflow Strategy Comparison

|approach|decision role|main strength|main risk|
|---|---|---|---|
|original_scripts|current end-to-end biological reference|uses all available information and currently supports downstream HMM/driver biology|stability may partially reflect repeated ST8-heavy covariance structure|
|InputQC_sensitivity|evidence layer for whether original input assumptions are safe|identifies whether depth, low-detection, and ST8 structure threaten interpretation|diagnostic rather than a complete alternative network decision by itself|
|balancednetwork|sampling-fair challenger and sensitivity control|directly tests whether modules survive after removing ST8 sample-count dominance|higher grey burden and smaller training set may discard real information|

## Remaining Gaps

|gap|status|next action|
|---|---|---|
|balanced downstream HMM rerun|missing; current HMM outputs are original/current pipeline|run HMM using balanced top3 module eigengenes in an isolated output directory|
|balanced driver/taxon importance rerun|missing; current importance outputs are original/current pipeline|rerun taxon importance and driver integration against balanced HMM/module outputs|
|balanced functional enrichment/state breakdown|missing for balanced; available for original/current pipeline|rerun functional/state summaries after balanced HMM|
|matched original-on-balanced-samples run|not yet run as a dedicated control|fit exp3 parameters on the balanced sample manifest and compare to balanced top3|
|module overlap original exp3 vs balanced top3|partly available through balanced_vs_original overlap, needs integrated summary|add per-module best-overlap table with taxon/function annotations|
|age-bin/core-stratified downstream validation|not yet summarized in method comparison|stratify key state/driver claims by age bin and core|

## Output Tables

- `tables/method_level_comparison.tsv`
- `tables/setting_level_comparison.tsv`
- `tables/module_distribution_comparison.tsv`
- `tables/module_stability_comparison.tsv`
- `tables/kme_topology_comparison.tsv`
- `tables/downstream_readiness_checklist.tsv`
- `tables/workflow_strategy_comparison.tsv`
- `tables/analysis_gap_matrix.tsv`
