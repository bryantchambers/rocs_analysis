# kME and Topology QC Strategy

Generated: 2026-05-13

Purpose: define the next NetworkQC step for comparing `baseline`, `opt5`, `exp4`, and `exp3` using WGCNA-native module membership and TOM topology metrics. This stays inside `networkQC` and does not change the main pipeline.

## Current Position

The preprocessing review found that read-depth structure remains visible across multiple input choices. That is a project-level caution, but it does not by itself choose a better WGCNA parameter set. For now, keep the current input fixed and ask a more local question:

Do the candidate module assignments form strong, coherent, internally connected modules?

The comparison set is:

|setting|power|deepSplit|mergeCutHeight|minModuleSize|current role|
|---|---:|---:|---:|---:|---|
|baseline|20|2|0.15|20|original/main-pipeline reference|
|opt5|12|1|0.20|20|best older 5-module-style option|
|exp4|12|3|0.20|20|runner-up, strongest preservation count|
|exp3|12|3|0.25|20|current best full-eval setting|

## What kME Adds

kME is the correlation between a taxon's abundance profile and a module eigengene. In plain terms, it asks: "does this taxon actually move with the module it was assigned to?"

This matters because WGCNA modules are not only hard clusters. A taxon can be assigned to a module but still be weakly attached to it. That is especially important here because downstream driver analysis depends on the biological meaning of module membership.

### kME Metrics

Compute these for each setting, separately for biological modules and technical modules:

- `median_assigned_kME`: typical strength of assigned module membership.
- `p05_assigned_kME`: lower-tail membership quality; low values mean some members barely belong.
- `frac_assigned_kME_lt_0_2`: fraction of weakly attached taxa.
- `frac_negative_assigned_kME`: fraction moving opposite to their assigned eigengene.
- `frac_assigned_is_max_kME`: fraction where assigned module is also the taxon's strongest module by kME.
- `median_kME_margin`: assigned kME minus next-best kME; low margins mean ambiguous membership.
- `n_strong_hubs`: taxa with assigned kME >= 0.70.
- `grey_max_bio_kME_median`: for grey taxa, median maximum kME to any biological module.
- `grey_rescuable_fraction`: grey taxa with max biological kME >= 0.50.

Interpretation:

- High assigned kME and high assigned-is-max are good.
- Low weak/negative membership is good.
- A high grey-rescuable fraction means the setting may be throwing biologically coherent taxa into grey.
- `exp3` should beat baseline here if the lower grey fraction is truly recovering coherent biological signal rather than just over-splitting.

## What Topology Adds

TOM asks whether two taxa share network neighbors, not just whether they correlate directly. It is the core WGCNA graph structure used for module detection.

Topology QC asks: "are within-module TOM connections stronger than between-module TOM connections?"

### TOM Metrics

Recompute TOM on the pooled training cores for each unique soft power:

- power 12 for `opt5`, `exp4`, and `exp3`.
- power 20 for `baseline`.

Then evaluate each setting's module labels on the corresponding TOM.

Compute:

- `within_tom_median`: median TOM for taxon pairs inside the same biological module.
- `between_tom_median`: median TOM for taxon pairs in different biological modules.
- `tom_separation_ratio`: within / between TOM median.
- `tom_silhouette_like`: `(within - between) / max(within, between)`.
- `within_edge_fraction_top_0_25pct`, `within_edge_fraction_top_0_5pct`, `within_edge_fraction_top_1pct`, `within_edge_fraction_top_2pct`.
- `modularity_top_0_5pct`: modularity of the thresholded TOM graph using WGCNA module labels.
- `largest_component_fraction_top_0_5pct`.
- `isolate_fraction_top_0_5pct`.

Interpretation:

- High within/between separation is good.
- High top-edge within-module fraction is good.
- High modularity is good, but only if not caused by nearly all taxa becoming isolates.
- Threshold sensitivity matters. A setting should not look good only at one arbitrary top-edge cutoff.

## Integration With Existing Full Evaluation

The final NetworkQC ranking should not replace the existing full-eval score. It should extend it.

Keep the existing full-eval score as the stability/preservation backbone because it already captures:

- bootstrap Jaccard stability,
- preservation in GeoB25202_R2,
- age-aligned eigengene concordance,
- core-balance sensitivity,
- grey burden and module count.

Add two new sub-scores:

- `kME_score`: module membership quality and low ambiguity.
- `topology_score`: TOM separation and robust threshold behavior.

Recommended final weighting:

- `full_eval_score`: 55%
- `kME_score`: 25%
- `topology_score`: 20%

Reasoning:

- Stability and preservation remain the main evidence.
- kME is the most WGCNA-native missing quality check.
- Topology is critical, but should not dominate because TOM threshold summaries can be sensitive to sparse network visualization choices.

## Decision Rules

Use these rules after the kME/topology run:

- Prefer `exp3` if it keeps the best overall score and does not show weak kME membership or poor TOM separation.
- Prefer `exp4` over `exp3` only if it materially improves kME/topology while staying close on grey burden and bootstrap/core-balance stability.
- Prefer `opt5` only if the expanded-module settings show clear over-fragmentation, weak kME margins, or poor TOM separation.
- Do not return to `baseline` unless all optimized settings fail membership/topology checks. Baseline's grey burden is too high to be competitive otherwise.

Suggested failure flags:

- biological median assigned kME < 0.50,
- biological p05 assigned kME < 0.10,
- biological negative assigned kME > 5%,
- assigned-is-max fraction < 70%,
- TOM separation ratio < 1.25,
- within-edge fraction at top 0.5% < 60%,
- grey-rescuable fraction > 25%.

These thresholds are review flags, not automatic rejection rules.

## Implementation Plan

Add one focused script:

`networkQC/scripts/10_kme_topology_review.R`

Inputs:

- `results/stage1/prokaryotes_vst.rds`
- `results/stage1/sample_metadata_stage1.tsv`
- `networkQC/results/full_eval/<setting_id>/module_assignments.tsv`
- `networkQC/results/full_eval/<setting_id>/module_eigengenes.tsv`
- `networkQC/results/tables/full_eval_setting_metric_summary.tsv`

Outputs:

- `networkQC/results/tables/kme_module_membership_summary.tsv`
- `networkQC/results/tables/kme_taxon_membership.tsv`
- `networkQC/results/tables/topology_quality_summary.tsv`
- `networkQC/results/tables/topology_threshold_sensitivity.tsv`
- `networkQC/results/tables/final_qc_integrated_ranking.tsv`
- `networkQC/results/figures/kme_topology_setting_heatmap.png`
- `networkQC/results/figures/kme_module_distributions_<setting>.png`
- `networkQC/results/figures/topology_threshold_sensitivity.png`
- `networkQC/results/KME_TOPOLOGY_QC_REPORT.md`

Runtime expectations:

- kME is cheap.
- TOM is the expensive part. Recompute TOM only twice: once for power 12 and once for power 20.
- Do not bootstrap this step initially. Bootstrap membership stability can be added later only if the first review is ambiguous.

## Final Decision Shape

The final report should conclude with one of three outcomes:

- `Keep exp3`: best combined stability, membership, and topology.
- `Switch to exp4`: similar stability but better membership/topology.
- `Hold opt5 as conservative`: expanded settings over-fragment or show weak membership.

The report should also carry forward the depth caveat from input evaluation:

All current candidates still show evidence of depth-related structure, so the selected network is the best available consensus network under the current input strategy, not proof that the depth issue has been fully solved.
