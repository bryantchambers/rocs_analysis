# Input Decision Report

- Generated: 2026-05-13 14:56:07 CEST
- Scope: focused decision review for `exp3` across four preprocessing inputs

## Plain-language summary

The current network uses a taxon-centered log transform. That means each taxon is compared to its own average across samples.
It is useful, but it is not the same thing as a sample-wise CLR, where each sample is centered against all taxa inside that sample.

Why this matters: if the transform leaves sample-wide read-depth structure in place, WGCNA may partly organize taxa around technical recovery rather than only around biology.

This review asks whether the current `exp3` modules still tell the same biological story when we switch the input frame.

## Decision table

|rank|input|class|score|mean overlap|grey %|PC1-depth|mean aligned ME corr|mean abs d18O shift|mean abs state shift|
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|
|1|Current exp3|current anchor|0.684|1.000|28.66|0.894|1.000|0.000|0.000|
|2|DESeq + length log|candidate alternative|0.452|0.633|37.62|0.877|0.897|0.033|0.088|
|3|Depth residualized|technical-control only|0.333|0.488|40.46|0.000|0.588|0.086|0.189|
|4|Sample CLR raw|robustness comparator|0.274|0.425|24.93|0.769|0.607|0.172|0.131|

## Interpretation

- Best overall usable input in this comparison: `current_taxon_centered_log`.
- `current exp3`: overlap 1.000, grey 28.66%, |PC1-depth| 0.894.
- `deseq_length_log`: overlap 0.633, grey 37.62%, |PC1-depth| 0.877.
- `sample_clr_raw`: overlap 0.425, grey 24.93%, |PC1-depth| 0.769.
- `log_depth_residualized`: overlap 0.488, grey 40.46%, |PC1-depth| 0.000.

Read these four inputs this way:

- `current exp3`: the continuity anchor. It preserves the existing network by definition and keeps the best overlap.
- `deseq_length_log`: the nearest non-current alternative. If it keeps similar module biology while cleaning technical structure, it is the first real replacement candidate.
- `sample_clr_raw`: the stronger compositional check. If it changes the network a lot, that does not prove it is wrong; it means the choice of compositional frame materially changes the network.
- `log_depth_residualized`: the control. If biology disappears only here, then depth is carrying part of the structure. This is useful diagnostically even if we never use it as the main input.

## Recommendation

Keep as the operational anchor for now, but the depth signal is still strong; treat `deseq_length_log` as the first fallback if downstream biology looks depth-driven. 

Do not treat `log_depth_residualized` as the production default unless a later review shows that the current and DESeq-based inputs are clearly dominated by technical structure.

## Outputs

- Summary: `networkQC/input_evaluation/results/decision/input_decision_summary.tsv`
- Module matching: `networkQC/input_evaluation/results/decision/module_match_table.tsv`
- Technical correlations: `networkQC/input_evaluation/results/decision/module_technical_correlations.tsv`
- Biological correlations: `networkQC/input_evaluation/results/decision/module_biological_correlations.tsv`
- Heatmap: `networkQC/input_evaluation/results/figures/input_decision_heatmap.png`
