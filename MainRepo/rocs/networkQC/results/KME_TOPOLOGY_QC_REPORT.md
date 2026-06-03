# kME and Topology QC Report

- Generated: 2026-05-13 15:22:44 CEST
- Scope: `baseline`, `opt5`, `exp4`, and `exp3`
- Progress log: `networkQC/results/kme_topology_review.log`

## Summary

This review adds two WGCNA-native checks to the existing full-eval ranking: kME module membership and TOM topology quality.

The input-depth artifact remains a standing caveat. This report chooses the best available module parameters under the current input strategy; it does not claim the depth issue is solved.

## Integrated Ranking

|rank|setting|full eval|kME|topology|integrated|grey %|bio median kME|TOM separation|flags|
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
|1|exp3|0.926|0.476|0.418|0.753|28.66|0.703|3.532|0|
|2|exp4|0.868|0.309|0.034|0.588|28.66|0.710|3.037|0|
|3|opt5|0.546|0.408|0.000|0.353|34.06|0.689|2.487|0|
|4|baseline|0.227|0.714|0.703|0.319|66.56|0.783|48.011|1|

## Recommendation

Outcome: `Keep exp3`.

Top setting: `exp3` with integrated score 0.753.

Decision rule: prefer `exp3` if it keeps the best combined score without weak kME or topology. Prefer `exp4` only if its kME/topology gain offsets its lower bootstrap/core-balance stability. Treat `opt5` as the conservative fallback if expanded settings look over-fragmented.

## Review Flags

Flags are warning signs, not automatic rejection rules. They mark settings that deserve inspection before being promoted into the main pipeline.

|setting|low median kME|low p05 kME|negative kME|low assigned-is-max|low TOM separation|low within edges|grey rescuable|total flags|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|exp3|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|0|
|exp4|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|0|
|opt5|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|0|
|baseline|FALSE|FALSE|FALSE|FALSE|FALSE|FALSE|TRUE|1|

## Outputs

- kME module summary: `networkQC/results/tables/kme_module_membership_summary.tsv`
- kME taxon table: `networkQC/results/tables/kme_taxon_membership.tsv`
- Topology summary: `networkQC/results/tables/topology_quality_summary.tsv`
- Threshold sensitivity: `networkQC/results/tables/topology_threshold_sensitivity.tsv`
- Integrated ranking: `networkQC/results/tables/final_qc_integrated_ranking.tsv`
- Heatmap: `networkQC/results/figures/kme_topology_setting_heatmap.png`
- Threshold plot: `networkQC/results/figures/topology_threshold_sensitivity.png`
