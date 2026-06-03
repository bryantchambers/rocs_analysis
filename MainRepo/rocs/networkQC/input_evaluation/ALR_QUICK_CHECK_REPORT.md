# ALR Quick Check Report

- Generated: 2026-05-13 12:26:26 CEST
- Scope: ALR-style input sensitivity for `exp3`, `exp4`, and `opt5`
- Preservation permutations: `50`

## Reference Choice

The preprint argues that ALR makes the reference explicit. Here, because we are building a network rather than testing a binary differential-abundance contrast, the reference was chosen as a technical-neutral denominator: high prevalence, nonzero median count, low association with total reads, age, and core.

- Single reference taxon: `S__GCA_003331785.1`
- Reference panel size: `20` taxa

## Summary

|setting|variant|grey_pct|mean_overlap|min_overlap|strong|moderate|pearson|spearman|rmse|PC1_depth|PC2_depth|drastic|
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
|exp3|alr_panel_reference|10.57|0.247|0.060|3|6|0.862|0.795|0.126|0.142|-0.781|yes|
|exp3|alr_single_reference|0.06|0.089|0.025|1|0|0.766|0.792|0.141|0.114|0.819|yes|
|exp4|alr_panel_reference|10.57|0.334|0.073|3|8|0.869|0.802|0.126|0.142|-0.781|yes|
|exp4|alr_single_reference|0.06|0.079|0.025|1|0|0.766|0.792|0.141|0.114|0.819|yes|
|opt5|alr_panel_reference|15.69|0.333|0.142|3|1|0.907|0.869|0.125|0.142|-0.781|yes|
|opt5|alr_single_reference|0.06|0.132|0.060|1|0|0.766|0.792|0.141|0.114|0.819|yes|

## Interpretation Guide

A good ALR result would reduce depth association across the leading PCs while keeping high overlap with the current best modules, high preservation, and good age-aligned concordance. If ALR strongly changes modules, it is informative but does not automatically solve the preprocessing problem; it means the reference frame controls the network geometry.

In this run, ALR reduced PC1-depth correlation but depth remained strong on PC2. The single-reference ALR also collapsed to nearly one non-grey module, and the panel-reference ALR changed module assignments substantially. This does not look like an immediate improvement over the previous input candidates.

## Output Map

- Reference candidates: `networkQC/input_evaluation/results/alr_quick_check/alr_reference_candidates.tsv`
- Summary: `networkQC/input_evaluation/results/alr_quick_check/alr_quick_check_summary.tsv`
- Module overlap: `networkQC/input_evaluation/results/alr_quick_check/alr_module_overlap_to_current.tsv`
- Preservation: `networkQC/input_evaluation/results/alr_quick_check/alr_preservation_all.tsv`
- Concordance: `networkQC/input_evaluation/results/alr_quick_check/alr_eigengene_concordance_all.tsv`
- Heatmap: `networkQC/input_evaluation/results/alr_quick_check/figures/alr_quick_check_heatmap.png`
