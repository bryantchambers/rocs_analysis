# NetworkQC Input Evaluation Rebuild Report

- Generated: 2026-05-13 11:27:26 CEST
- Scope: focused preprocessing/depth sensitivity for `exp3`, `exp4`, and `opt5`
- Preservation permutations per sensitivity run: `100`

## Input Variants

- `current_taxon_centered_log`: current main-pipeline WGCNA input, taxon-centered log raw counts.
- `sample_clr_raw`: sample-wise CLR on raw filtered counts with pseudocount 0.5.
- `deseq_length_log`: DESeq2 poscounts + reference-length normalized log counts, taxon-centered.
- `log_depth_residualized`: current log matrix after removing linear log-total-read signal per taxon.

## Summary Table

|setting|variant|grey_pct|mean_overlap|min_overlap|strong|moderate|pearson|spearman|rmse|PC1_depth|drastic|
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
|exp3|current_taxon_centered_log|28.66|1.000|1.000|6|2|0.910|0.726|0.078|-0.894|no|
|exp3|deseq_length_log|37.62|0.633|0.528|6|4|0.893|0.715|0.087|-0.877|no|
|exp3|log_depth_residualized|40.46|0.488|0.089|3|14|0.787|0.778|0.138|-0.000|yes|
|exp3|sample_clr_raw|24.93|0.425|0.174|7|11|0.867|0.831|0.137|0.769|yes|
|exp4|current_taxon_centered_log|28.66|1.000|1.000|7|2|0.907|0.739|0.085|-0.894|no|
|exp4|deseq_length_log|37.62|0.508|0.183|7|7|0.880|0.746|0.101|-0.877|yes|
|exp4|log_depth_residualized|40.46|0.458|0.145|3|14|0.787|0.778|0.138|-0.000|yes|
|exp4|sample_clr_raw|24.93|0.449|0.174|7|11|0.867|0.831|0.137|0.769|yes|
|opt5|current_taxon_centered_log|34.06|1.000|1.000|5|0|0.912|0.726|0.092|-0.894|no|
|opt5|deseq_length_log|43.02|0.534|0.173|6|1|0.892|0.762|0.100|-0.877|yes|
|opt5|log_depth_residualized|46.74|0.496|0.169|3|7|0.746|0.769|0.141|-0.000|yes|
|opt5|sample_clr_raw|39.51|0.497|0.302|6|2|0.907|0.859|0.133|0.769|yes|

## Initial Interpretation

At least one non-current variant crossed the pre-set drastic-change thresholds. Review these before locking the current network:

- `exp3 / log_depth_residualized`: mean overlap 0.488, min overlap 0.089, grey delta 11.80, strong biological modules 3, Pearson 0.787
- `exp3 / sample_clr_raw`: mean overlap 0.425, min overlap 0.174, grey delta -3.73, strong biological modules 7, Pearson 0.867
- `exp4 / deseq_length_log`: mean overlap 0.508, min overlap 0.183, grey delta 8.96, strong biological modules 7, Pearson 0.880
- `exp4 / log_depth_residualized`: mean overlap 0.458, min overlap 0.145, grey delta 11.80, strong biological modules 3, Pearson 0.787
- `exp4 / sample_clr_raw`: mean overlap 0.449, min overlap 0.174, grey delta -3.73, strong biological modules 7, Pearson 0.867
- `opt5 / deseq_length_log`: mean overlap 0.534, min overlap 0.173, grey delta 8.96, strong biological modules 6, Pearson 0.892
- `opt5 / log_depth_residualized`: mean overlap 0.496, min overlap 0.169, grey delta 12.69, strong biological modules 3, Pearson 0.746
- `opt5 / sample_clr_raw`: mean overlap 0.497, min overlap 0.302, grey delta 5.45, strong biological modules 6, Pearson 0.907

For `exp3`, compare the non-current rows against `current_taxon_centered_log`. If module overlap stays high and preservation/concordance remain strong, the current `exp3` decision is robust to input preprocessing. If overlap collapses or depth-controlled variants change the module structure, rerun the broader parameter sweep on the preferred corrected input.

## Output Map

- Input variants: `networkQC/input_evaluation/results/inputs/*.rds`
- Summary table: `networkQC/input_evaluation/results/sensitivity/input_sensitivity_summary.tsv`
- Module overlap details: `networkQC/input_evaluation/results/sensitivity/module_overlap_to_current.tsv`
- Preservation details: `networkQC/input_evaluation/results/sensitivity/preservation_all.tsv`
- Concordance details: `networkQC/input_evaluation/results/sensitivity/eigengene_concordance_all.tsv`
- Heatmap: `networkQC/input_evaluation/results/figures/input_sensitivity_heatmap.png`
