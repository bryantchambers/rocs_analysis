# Network QC Decision Report

- Generated: 2026-05-07 10:18:57 CEST
- Inputs: parameter sweep + cross-core eigengene concordance

## Scoring Strategy

The decision score is weighted toward lower grey burden and realistic module count,
with module-size sanity as a secondary factor.

- `score_grey` (45%): lower `grey_pct` is better
- `score_module_count` (35%): closer to 5 non-grey modules is better
- `score_module_size` (20%): larger median non-grey module size is better

## Cross-core Concordance Benchmark

- Mean Pearson r: 0.526
- Mean Spearman rho: 0.483
- Mean RMSE: 0.091

## Top 3 Parameter Sets

|rank|power|deepSplit|mergeCutHeight|minModuleSize|non_grey_modules|grey_pct|module_size_median|decision_score|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|1|12|1|0.10|20|5|34.06|259.0|0.925|
|2|12|1|0.10|30|5|34.06|259.0|0.925|
|3|12|1|0.15|20|5|34.06|259.0|0.925|

## Recommendation

Run full stability diagnostics (`02b`-style bootstrap) on these top 3 settings,
then pick the one with best joint behavior across:
1) bootstrap Jaccard stability, 2) biological-module preservation,
3) age-aligned eigengene concordance, and 4) downstream biological consistency.
