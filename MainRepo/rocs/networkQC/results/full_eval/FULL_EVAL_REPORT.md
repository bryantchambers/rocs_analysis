# NetworkQC Full Evaluation Report

- Generated: 2026-05-13 09:54:24 CEST
- Scope: baseline + top optimized sweep settings + expansion shortlist settings (full preservation + bootstrap)

## Ranking Table

|rank|setting_id|power|deepSplit|mergeCutHeight|minModuleSize|grey_pct|mean_bootstrap_jaccard|bio_pres_strong|mean_concordance_pearson|mean_balanced_jaccard|final_score|
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|1|exp3|12|3|0.25|20|28.66|0.400|6|0.857|0.517|0.926|
|2|exp4|12|3|0.20|20|28.66|0.378|7|0.859|0.459|0.868|
|3|opt5|12|1|0.20|20|34.06|0.395|5|0.714|0.456|0.546|
|4|exp1|12|2|0.25|30|33.22|0.390|4|0.741|0.426|0.540|
|5|opt3|12|1|0.15|20|34.06|0.392|5|0.714|0.447|0.529|
|6|opt4|12|1|0.15|30|34.06|0.378|5|0.714|0.511|0.509|
|7|opt2|12|1|0.10|30|34.06|0.378|5|0.714|0.511|0.508|
|8|opt1|12|1|0.10|20|34.06|0.386|5|0.714|0.447|0.505|
|9|exp2|12|2|0.25|20|29.10|0.338|5|0.739|0.463|0.433|
|10|baseline|20|2|0.15|20|66.56|0.351|4|0.847|0.349|0.227|
|11|exp5|16|1|0.15|20|48.14|0.330|5|0.757|0.402|0.209|

## Best Setting Recommendation

- Best setting: `exp3`
- Parameters: power=12, deepSplit=3, mergeCutHeight=0.25, minModuleSize=20
- Grey fraction: 28.66%
- Bootstrap stability (mean module median Jaccard): 0.400
- Biological preservation strong modules: 6
- Core-balance robustness (mean balanced Jaccard): 0.517

## Expanded Stability Review

Recommendation: use `exp3` as the current best WGCNA parameter set.

- `exp3` parameters: power=12, deepSplit=3, mergeCutHeight=0.25, minModuleSize=20
- It has the lowest grey burden among the top settings (28.66%), 8 non-grey modules, 6 strong and 2 moderate biological modules, mean bootstrap Jaccard 0.400, and mean balanced Jaccard 0.517.
- Biological-module medians for `exp3`: Jaccard median 0.406, Jaccard p05 0.139, Jaccard p95 0.585, Zsummary 11.70, Zdensity 15.70, Zconnectivity 8.12, Pearson 0.913, Spearman 0.770, RMSE 0.066.
- Runner-up `exp4` has slightly stronger preservation count (7 strong, 2 moderate) and similar correlation, but lower mean bootstrap Jaccard 0.378 and lower mean balanced Jaccard 0.459.
- The previous best 5-module option `opt5` remains reasonable but retains more grey taxa (34.06%) and has lower age-aligned Pearson concordance 0.714.
- The original baseline is not competitive here: grey taxa remain 66.56% and mean balanced Jaccard is 0.349.

## Output Map

- Per-setting outputs: `networkQC/results/full_eval/<setting_id>/`
- Queued settings: `networkQC/results/full_eval/settings_to_evaluate.tsv`
- Ranked summary: `networkQC/results/full_eval/all_settings_ranked.tsv`
- Graph plots by setting: `networkQC/results/figures/full_graph_<setting_id>.png`
- Leiden graph plot: `networkQC/results/figures/full_graph_leiden.png`
- Setting-level heatmap: `networkQC/results/figures/full_eval_setting_metric_heatmap.png`
- Module-level heatmap: `networkQC/results/figures/full_eval_module_metric_heatmap.png`
- Setting-level metric table: `networkQC/results/tables/full_eval_setting_metric_summary.tsv`
- Module-level metric table: `networkQC/results/tables/full_eval_module_metric_summary.tsv`

## Interpretation Notes

1. `mean_bootstrap_jaccard` captures module reproducibility under resampling.
2. `bio_pres_strong` prioritizes non-technical module preservation in R1->R2.
3. `mean_balanced_jaccard` estimates how sensitive module formation is to core-size imbalance.
4. Choose settings with strong performance across all three, not just low grey.
