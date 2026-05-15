# InputQC Research Plan: Sequencing-Depth Signal and Biological Recovery

Generated: 2026-05-15

Purpose: create a contained plan for investigating why the leading axis of the WGCNA input is strongly associated with sequencing depth, and identify defensible input transformations / correction strategies before changing the main pipeline.

## Current Starting Point

We found that the main WGCNA input is not a true sample-wise CLR. In `scripts/01_data_prep.R`, the transform is:

```r
clr_mat <- log(count_mat + 0.5)
clr_mat <- clr_mat - rowMeans(clr_mat)
```

Because `count_mat` is taxa x samples, this centers each taxon across samples. In plain language, each taxon is put on an "above or below its own average" scale. That can be useful for expression-style correlation, but it does not center each sample by its own total composition.

In NetworkQC, this mattered because the current input's PC1 was strongly correlated with total reads:

- current taxon-centered log input: `PC1 ~ log_total_reads r = -0.894`
- `deseq_length_log`: `r = -0.877`
- `sample_clr_raw`: `r = 0.769`
- `log_depth_residualized`: `r ~ 0`, but module structure changed substantially

Interpretation: sequencing depth / DNA preservation is not just a superficial normalization issue. It may be a major axis of the data itself.

## InputQC Workspace

This folder is a contained subanalysis and should not write into the main pipeline outputs unless explicitly promoted later.

Copied inputs:

- `InputQC/source_copies/01_data_prep.main_pipeline.R`
- `InputQC/scripts/01_data_prep_input_variants.R`
- `InputQC/scripts/02_eval_input_sensitivity.R`
- `InputQC/scripts/03_alr_quick_check.R`
- existing copied reports and results from `networkQC/input_evaluation`

The first implementation goal is not to pick a final normalization. It is to understand the failure modes well enough to choose the next serious test.

## Research Questions

1. Is depth acting as a purely technical artifact, or is it confounded with real ancient DNA preservation / age / core / sediment chemistry?
2. Does the depth signal remain after correct sample-wise compositional transforms?
3. Does the signal reflect total reads, damage-authenticated reads, reference length, library concentration, fragment length, core identity, or age?
4. Can we reduce depth structure without destroying climate/state/module biology?
5. Which correction is appropriate for a downstream correlation/network input rather than differential abundance testing?

## Literature Anchors and Practical Meaning

WGCNA guidance:

- WGCNA is correlation-based, so sample-wise scaling and technical shifts matter.
- The WGCNA FAQ recommends checking sample quantiles and correcting systematic sample shifts when they are technical rather than biological.
- WGCNA itself does not solve normalization problems; the input matrix must already be biologically meaningful.
- Source: WGCNA FAQ, https://edo98811.github.io/WGCNA_official_documentation/faq.html
- Source: Langfelder and Horvath 2008, https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559

Compositional microbiome guidance:

- Microbiome/metagenomic count tables are compositional because sample totals constrain observed counts.
- CLR/log-ratio approaches are common, but zeros, pseudocounts, sparse taxa, and sample heterogeneity can still create problems.
- Network inference on microbiome data should be interpreted as association structure, not direct ecological interaction.
- Source: Gloor et al. 2017, https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
- Source: NetCoMi, https://pmc.ncbi.nlm.nih.gov/articles/PMC8293835/
- Source: SPIEC-EASI, https://pmc.ncbi.nlm.nih.gov/articles/PMC4423992/

Variance stabilization:

- DESeq2 VST transforms normalized counts to reduce the mean-variance relationship and is useful for clustering / visualization / machine learning inputs.
- VST normalizes with respect to size factors or normalization factors, so it is a serious candidate for network input testing.
- It is not automatically compositional, and DESeq2 assumptions may not perfectly match ancient metagenomic data.
- Source: DESeq2 vignette, https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

Sampling fraction / bias correction:

- ANCOM-BC estimates sampling-fraction bias and is designed around microbiome compositional bias.
- It is primarily a differential-abundance method, not a direct network-input transform, but its framing is useful: library size is not necessarily the same as sampling fraction.
- Source: Lin and Peddada 2020, https://www.nature.com/articles/s41467-020-17041-7
- Source: ANCOM-BC2 vignette, https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

Unwanted variation / batch correction:

- Known technical variables can be adjusted directly; unknown technical variation can be estimated with SVA/RUV-like methods.
- For this project, these methods must be used carefully because depth, age, core, and preservation may be partly biological and partly technical.
- Over-correction could remove real climate/preservation biology.
- Source: sva Bioconductor package, https://bioconductor.org/packages/sva/
- Source: RUVSeq vignette, https://new.bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html

Ancient DNA caution:

- Ancient DNA has true preservation/damage structure: fragmentation, post-mortem damage, contamination, and variable recovery.
- Therefore, "depth correction" is not automatically correct. Some depth-associated variation may reflect preservation context, sediment conditions, or authentic historical signal.
- Source: Nature Reviews Genetics ancient genomics review, https://www.nature.com/articles/nrg3935
- Source: BMC Genomics contamination/aDNA review context, https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-2020-07229-y

## Candidate Strategy Families

### 1. Correct sample-wise CLR and pseudocount sensitivity

Question: does a correct sample-wise CLR reduce depth effects without breaking biological modules?

Tests:

- sample-wise CLR using pseudocounts `0.5`, `1`, and possibly multiplicative zero replacement.
- CLR after prevalence/abundance filters of increasing strictness.
- CLR on proportions rather than raw counts, only as a sanity check.
- robust PCA / ordinary PCA correlation with depth, age, core, fragment length, and library concentration.

Risks:

- CLR can still carry depth effects if zeros and detection limits are depth-dependent.
- Pseudocount choice can dominate rare taxa.

### 2. Variance-stabilized count models

Question: does a count-model transform reduce mean/depth variance better than log/CLR?

Tests:

- DESeq2 `vst()` using `poscounts` size factors.
- DESeq2 VST with reference-length normalization factors already used for EMP.
- rlog only if runtime is acceptable; likely slower and not needed first.
- compare PCA-depth, mean-SD trend, sample distances, and module stability.

Risks:

- DESeq2 assumes most features are not systematically changing in ways that violate size-factor assumptions.
- VST may reduce variance structure but not solve compositional bias.

### 3. Microbiome-specific size-factor alternatives

Question: do microbiome-specific normalization factors reduce depth / sampling-fraction structure better than DESeq2?

Potential methods:

- GMPR-style size factors for zero-inflated microbiome data.
- CSS / metagenomeSeq-like normalization.
- TMM / upper-quartile as sensitivity checks.

Tests:

- transform normalized counts with log1p or VST-like scaling.
- compare PCA-depth, depth residuals, module stability, and biological signal retention.

Risks:

- Some methods are designed for differential abundance, not correlation networks.
- Normalization can improve sample comparability while still leaving compositional artifacts.

### 4. Known technical covariate regression

Question: can we remove technical axes while preserving age/climate/core biology?

Candidate covariates:

- log total reads
- initial reads
- dereplicated reads
- library concentration
- average initial fragment length
- average dereplicated fragment length
- core
- flowcell if usable
- damage metrics if available in upstream damage table

Tests:

- regress only clearly technical variables from transformed taxon matrix.
- compare three designs:
  - depth-only residualization
  - depth + library/fragment metrics
  - depth + technical metrics while preserving climate/age terms
- evaluate whether climate/state associations disappear.

Risks:

- Depth and preservation may be biologically entangled with age/core.
- Regressing depth can remove real preservation-linked ecological structure.
- A residualized matrix should initially be treated as a control, not the main input.

### 5. Latent unwanted-variation estimation

Question: are there hidden technical factors beyond measured depth?

Potential methods:

- SVA / svaseq on transformed counts.
- RUV-style unwanted factor estimation using negative-control taxa if we can define them.
- PCA/PEER-style latent factor removal as exploratory sensitivity.

Tests:

- identify latent factors and correlate them with depth, core, age, fragment length, library concentration, and climate proxies.
- remove 1 factor at a time and evaluate biology retention.

Risks:

- If latent factors align with real biology, removing them will erase the signal we care about.
- We need negative controls or strong assumptions before using RUV-style correction as a main input.

### 6. Filtering and sample inclusion sensitivity

Question: is the depth axis driven by low-quality samples or sparse taxa?

Tests:

- sample-depth thresholds: remove very low total-read samples.
- damage-authenticated-read thresholds.
- prevalence thresholds: 10, 15, 20 samples.
- abundance thresholds: remove ultra-rare taxa even if prevalent.
- core-specific checks: confirm no single core drives the depth axis.

Risks:

- Filtering can improve stability but reduce the temporal record.
- Removing low-depth samples may bias the age/core structure.

### 7. Biological-signal retention checks

Every candidate input must be judged by both artifact reduction and biology retention.

Required outputs for each candidate:

- PC1/PC2/PC3 correlations with total reads and other technical variables.
- PC correlations with age, core, d18O, SST.
- sample ordination plots colored by depth, age, core, d18O, SST.
- module stability if used for WGCNA.
- HMM-state or climate eigengene concordance if modules are rebuilt.
- preservation and kME/topology if a candidate becomes serious.

Decision principle:

Do not select the input with the lowest depth correlation if it destroys the biological structure. Select the input that best reduces technical dominance while preserving stable, interpretable climate/state/module biology.

## Proposed Stepwise Workflow

### Phase 1: Diagnostic inventory

Build `InputQC/scripts/00_input_diagnostics.R`.

Outputs:

- `InputQC/results/tables/sample_technical_covariates.tsv`
- `InputQC/results/tables/technical_covariate_correlations.tsv`
- `InputQC/results/tables/input_pc_associations.tsv`
- `InputQC/results/figures/input_pc_technical_heatmap.png`
- `InputQC/results/figures/sample_ordination_by_depth_core_age.png`
- `InputQC/INPUT_QC_DIAGNOSTIC_REPORT.md`

Goal:

Map the technical/biological confounding structure before applying more corrections.

### Phase 2: Candidate input generation

Build `InputQC/scripts/01_generate_candidate_inputs.R`.

Candidate inputs:

- current taxon-centered log
- sample-wise CLR variants
- DESeq2 VST with poscounts
- DESeq2 VST + reference-length normalization
- GMPR/CSS/TMM/log variants if packages are available
- depth-only residualized control
- technical-covariate residualized control
- one or two filtered-data variants

Goal:

Create a controlled set of input matrices with consistent sample/taxon IDs.

### Phase 3: Input scoring without WGCNA

Build `InputQC/scripts/02_score_candidate_inputs.R`.

Score each input on:

- depth association reduction
- technical covariate reduction
- age/climate/core interpretability
- mean-variance behavior
- sample distance stability
- taxa retention

Goal:

Avoid rerunning full WGCNA for every speculative transform.

### Phase 4: WGCNA mini-evaluation on shortlisted inputs

For the top 2-3 candidate inputs only:

- rebuild WGCNA with `exp3` parameters.
- compare module overlap to current `exp3`.
- compute preservation, eigengene concordance, kME, and topology.
- treat residualized inputs as controls unless clearly superior.

Goal:

Only promote an input if it improves technical behavior without breaking the network.

### Phase 5: Decision conference

Before changing the main pipeline:

- review the diagnostic report.
- decide whether to keep current input with caveat, switch to a transformed input, or rerun broader NetworkQC on a new input.
- update `CODEX.md`, `BryantsNotes.md`, and pipeline docs.

## Initial Hypotheses

1. The depth signal will not disappear completely because sequencing depth is likely tied to ancient DNA preservation.
2. Sample-wise CLR will reduce some compositional artifacts but may not solve detection-limit effects from sparse/low-depth samples.
3. DESeq2 VST with reference-length normalization may be the most practical next serious candidate because it addresses size factors, reference lengths, and mean-variance behavior.
4. Direct depth residualization will probably remove too much structure and should remain a control unless biological signals are robust afterward.
5. Filtering low-quality samples/taxa may matter as much as the transformation.

## Questions for Bryant / Project Decision

1. Are we willing to remove low-depth samples if they drive PC1, even if this reduces age/core coverage?
2. Which biological signals are non-negotiable to preserve during input correction: HMM states, d18O/SST relationships, module preservation in R2, or functional EMP/TEA links?
3. Are there upstream damage-authentication metrics beyond `is_dmg`, fragment length, and read counts that we should add to the technical covariate table?
4. Should `core` be treated as biology, technical structure, or both? This matters for residualization.
5. Should we prioritize a conservative input that preserves the current `exp3` story, or a stricter input that reduces depth even if modules change?

## Sources

- WGCNA FAQ: https://edo98811.github.io/WGCNA_official_documentation/faq.html
- WGCNA original paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559
- Gloor et al. compositional microbiome data: https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
- NetCoMi microbiome network review/toolkit: https://pmc.ncbi.nlm.nih.gov/articles/PMC8293835/
- SPIEC-EASI: https://pmc.ncbi.nlm.nih.gov/articles/PMC4423992/
- DESeq2 vignette: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- ALDEx2 vignette: https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html
- ANCOM-BC paper: https://www.nature.com/articles/s41467-020-17041-7
- ANCOM-BC2 vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
- sva Bioconductor package: https://bioconductor.org/packages/sva/
- RUVSeq vignette: https://new.bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html
- Ancient genomics review: https://www.nature.com/articles/nrg3935
