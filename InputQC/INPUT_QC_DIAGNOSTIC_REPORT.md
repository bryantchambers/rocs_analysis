# InputQC Diagnostic Report

- Generated: 2026-05-15 09:35:11 CEST
- Scope: diagnostic inventory before new normalization/correction choices
- Progress log: `InputQC/results/input_qc_diagnostics.log`

## Main finding

The current input still has a strong PC1 association with log total reads: `r = -0.894`.

This confirms that the next step should be diagnosis before correction. In ancient DNA, depth can be both technical and preservation-linked, so removing it blindly could erase real structure.

## Top PC associations

|variant|PC|covariate|class|Pearson r|Spearman rho|
|---|---|---|---|---:|---:|
|current_taxon_centered_log|PC1|detected_taxa|technical|-0.987|-0.980|
|deseq_length_log|PC1|detected_taxa|technical|-0.911|-0.917|
|current_taxon_centered_log|PC1|log_total_reads|technical|-0.894|-0.906|
|deseq_length_log|PC1|log_total_reads|technical|-0.877|-0.792|
|deseq_length_log|PC1|shannon_raw|technical|-0.855|-0.795|
|sample_clr_raw|PC1|detected_taxa|technical|0.844|0.815|
|current_taxon_centered_log|PC1|log_derep|technical|-0.818|-0.822|
|deseq_length_log|PC1|log_derep|technical|-0.816|-0.716|
|current_taxon_centered_log|PC1|shannon_raw|technical|-0.776|-0.774|
|sample_clr_raw|PC1|log_total_reads|technical|0.769|0.703|
|sample_clr_raw|PC1|shannon_raw|technical|0.767|0.771|
|current_taxon_centered_log|PC1|total_reads|technical|-0.740|-0.906|

## Method assumption register

The goal is not to reject a method because it was designed elsewhere. The goal is to track what assumptions it may impose on ancient DNA data.

|method family|useful for|key assumption risk|current stance|
|---|---|---|---|
|Taxon-centered log current input|Continuity with existing WGCNA and comparable module story|Does not remove sample-wide depth effects; can preserve technical axes|Anchor / current reference|
|Sample-wise CLR|Compositional sample centering and log-ratio geometry|Zeros and pseudocounts can dominate rare taxa; detection limits remain depth-dependent|Required stress test, not automatically better|
|DESeq2 VST / rlog|Mean-variance stabilization and size-factor normalization|Built for RNA-seq-like negative-binomial count behavior; assumes size factors are meaningful for this system|Serious next candidate, but assumptions must be tracked|
|Microbiome size factors: GMPR/CSS/TMM|Zero-inflated microbiome-like count normalization|Often designed for microbiome differential abundance, not necessarily correlation-network inputs|Worth testing if packages are available; not first-line without diagnostics|
|Known covariate residualization|Testing whether measured technical factors explain PCs|May erase real ancient-DNA preservation biology if depth/fragment metrics are biologically entangled|Control/sensitivity before production use|
|SVA/RUV latent factor removal|Testing unknown unwanted variation|Latent factors may represent true climate/core biology rather than nuisance variation|Exploratory only until negative controls or clear nuisance factors are defined|
|Filtering low-depth samples/taxa|Testing whether low-quality observations drive PC1|Can improve quality but may bias age/core coverage and remove informative ancient samples|Likely important; must report sample/age/core loss|

## Recommended next move

Use these diagnostics to choose a small set of candidate corrections. The first serious candidates should be:

- DESeq2 VST with poscounts and reference-length normalization, as a count-model variance-stabilized input.
- sample-wise CLR with stricter filtering and pseudocount sensitivity, as the compositional baseline.
- depth/technical residualization as a control only, not as the default production input.
- sample/taxon filtering sensitivity, because low-depth samples and sparse taxa may drive PC1 as much as the transform does.

Do not promote any corrected input unless it reduces technical dominance while preserving climate/state/module biology.

## Outputs

- Covariates: `InputQC/results/tables/sample_technical_covariates.tsv`
- Covariate correlations: `InputQC/results/tables/technical_covariate_correlations.tsv`
- PC associations: `InputQC/results/tables/input_pc_associations.tsv`
- Top PC associations: `InputQC/results/tables/input_pc_top_associations.tsv`
- Method assumptions: `InputQC/results/tables/method_assumption_register.tsv`
- Heatmap: `InputQC/results/figures/input_pc_technical_heatmap.png`
- Ordination by depth: `InputQC/results/figures/ordination_current_by_depth.png`
- Ordination by core: `InputQC/results/figures/ordination_current_by_core.png`
