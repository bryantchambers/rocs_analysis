# Network QC Report

- Generated: 2026-05-06 17:01:29 CEST
- Scope: WGCNA QC baseline, parameter sweep, Leiden comparison

## 1. Baseline WGCNA QC

Key metrics:

|metric|value|
|---|---:|
|n_samples_train|189|
|n_taxa|1797|
|n_modules_non_grey|5|
|grey_pct|66.5553700612131|
|top_edge_threshold|0.107868624603524|

## 2. Parameter Sweep Summary

- Total combinations tested: 72
- Successful fits: 72
- Lowest grey_pct candidate: power=12, deepSplit=3, mergeCutHeight=0.10, minModuleSize=20, grey_pct=28.66

## 3. Cross-core Eigengene Concordance

|core_a|core_b|pearson_mean|spearman_mean|rmse_mean|
|---|---|---:|---:|---:|
|GeoB25202_R1|GeoB25202_R2|0.847|0.696|0.097|
|ST13|GeoB25202_R1|0.807|0.718|0.047|
|ST13|GeoB25202_R2|0.757|0.705|0.102|
|ST8|GeoB25202_R1|0.290|0.279|0.085|
|ST8|GeoB25202_R2|0.272|0.284|0.114|
|ST8|ST13|0.183|0.217|0.101|

## 4. Leiden Comparison

Resolution scan:

|resolution|n_modules|module_size_median|largest_module|
|---:|---:|---:|---:|
|0.05|267|1.0|929|
|0.10|269|1.0|742|
|0.20|269|1.0|726|
|0.40|270|1.0|707|

Selected Leiden run:

|metric|value|
|---|---|
|selected_power|12|
|tom_threshold|0.05|
|selected_resolution|0.05|
|leiden_modules|267|

WGCNA-to-Leiden dominant mapping purity:

|module_wgcna|module_leiden|N|purity|
|---|---|---:|---:|
|grey|L1|545|1.000|
|turquoise|L1|237|1.000|
|blue|L3|133|1.000|
|brown|L3|84|1.000|
|yellow|L1|76|1.000|
|green|L1|65|1.000|

## 5. Next Actions

1. Use sweep results to choose 2-3 candidate WGCNA settings with lower grey burden.
2. Re-run stability (`02b`) for each candidate and compare Jaccard + concordance.
3. Compare downstream outcomes (07b/09/10/11/12) under WGCNA-best vs Leiden-best module sets.
