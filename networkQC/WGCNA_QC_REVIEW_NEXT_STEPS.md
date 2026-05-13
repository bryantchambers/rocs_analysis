# WGCNA / NetworkQC Review and Next-Step Notes

Generated: 2026-05-13

Purpose: capture the current QC interpretation, literature-grounded concerns, and reasonable next checks before we lock WGCNA module parameters for downstream biological analysis.

## Current Decision State

The expanded NetworkQC run supports `exp3` as the best current WGCNA parameter set:

- `power = 12`
- `deepSplit = 3`
- `mergeCutHeight = 0.25`
- `minModuleSize = 20`

Why `exp3` currently wins:

- Lowest grey burden among top settings: `28.66%`
- Good module granularity: `8` non-grey modules
- Biological preservation: `6` strong + `2` moderate biological modules
- Best mean bootstrap Jaccard in the expanded comparison: `0.400`
- Best mean balanced Jaccard: `0.517`
- Strong age-aligned eigengene concordance: Pearson `0.857`, Spearman `0.740`, RMSE `0.089`

`exp4` is the main runner-up. It has one more strongly preserved biological module, but weaker bootstrap/core-balance stability. I would only prefer `exp4` if we decide preservation count matters more than resampling stability.

## What We Already Do Well

The QC work now covers several important axes:

- Consensus WGCNA across training cores: ST8, ST13, GeoB25202_R1.
- GeoB25202_R2 held out for validation/preservation.
- Preservation testing with final-level permutations.
- Bootstrap module overlap stability using best-matched Jaccard.
- Age-aligned eigengene concordance using Pearson, Spearman, and RMSE.
- Core-balance sensitivity by downsampling training cores to equal sample counts.
- Parameter sweep beyond the original 5-module prior.
- Heatmaps summarizing setting-level and module-level stability/preservation/correlation metrics.
- Graph layout diagnostics showing that the earlier "shotgun" graph view was mostly an isolate-plotting artifact at the top 0.5% TOM edge threshold.

This is a defensible stability framework. The remaining question is not "did we do enough QC at all?" It is whether we have missed a small number of checks that could materially change confidence in downstream biology.

## Main Concern Found During Review

The largest issue is preprocessing-related.

In `scripts/01_data_prep.R`, the matrix saved as `prokaryotes_vst.rds` is described as CLR, but the operation is:

```r
clr_mat <- log(count_mat + 0.5)
clr_mat <- clr_mat - rowMeans(clr_mat)
```

Because `count_mat` is taxa x samples, `rowMeans(clr_mat)` centers each taxon across samples. That is useful for expression-style correlation, but it is not standard sample-wise CLR. A standard sample-wise CLR would subtract the sample-wise mean log abundance across taxa.

Practical consequence:

- WGCNA is currently close to being run on log raw counts with a pseudocount, then taxon-centered.
- Pearson correlation is invariant to taxon-wise centering/scaling, so this does not remove sample-wise library/depth effects.
- I checked the current WGCNA matrix and found PC1 strongly correlated with log total reads: about `r = -0.894`.
- Total reads vary from about `4,586` to `3,774,170`.
- Core medians also differ: ST8 median total reads are lower than GeoB/ST13.

This does not automatically invalidate the current network, but it is a genuine safety concern. Before locking the modules, we should test whether `exp3` remains stable under a more compositionally appropriate or depth-controlled input.

## Literature Anchors

WGCNA core guidance:

- WGCNA is built around correlation networks, module eigengenes, module membership/kME, signed networks, TOM, and module preservation.
- WGCNA documentation recommends signed networks in many settings, robust correlation (`bicor`) when outliers are plausible, and checking sample normalization/batch/technical drivers.
- The WGCNA FAQ explicitly notes that sample-wise scaling factors matter for correlation networks.

Useful sources:

- WGCNA manual: https://cran.r-universe.dev/WGCNA/doc/manual.html
- WGCNA FAQ: https://edo98811.github.io/WGCNA_official_documentation/faq.html
- Langfelder and Horvath 2008: https://doi.org/10.1186/1471-2105-9-559

Module preservation:

- Langfelder et al. 2011 emphasize that overlap/cross-tabulation alone is insufficient because it can miss preserved connectivity patterns.
- Preservation should evaluate density and connectivity, which is what `modulePreservation()` is designed to do.
- Standard rough interpretation: `Zsummary < 2` weak/no evidence, `2-10` moderate evidence, `>10` strong evidence.

Useful source:

- Langfelder et al. 2011: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001057

Microbiome/metagenomic network caution:

- Microbiome/metagenomic count data are compositional; direct correlation on unadjusted counts can induce artifacts.
- CLR/log-ratio approaches are common, but correlation of transformed compositional data still needs caution.
- Co-occurrence networks should be interpreted as association networks, not direct ecological interaction maps.
- Network methods such as SparCC, SPIEC-EASI, CCLasso, and NetCoMi exist because microbiome network inference has compositional and sparsity challenges.

Useful sources:

- Gloor et al. 2017: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full
- NetCoMi: https://academic.oup.com/bib/article/22/4/bbaa290/6017455
- Hirano and Takemoto 2019: https://link.springer.com/article/10.1186/s12859-019-2915-1

## What Our Current Metrics Measure

Overlap stability:

- Bootstrap Jaccard asks whether module membership is recovered when samples are resampled.
- Strength: directly tests assignment reproducibility.
- Weakness: can penalize split/merge behavior even when the underlying topology is similar.

Preservation:

- `Zsummary`, `Zdensity`, and `Zconnectivity` ask whether modules are preserved in GeoB25202_R2.
- Strength: evaluates connectivity structure, not just module overlap.
- Weakness: preservation is only R1 -> R2 right now, not all possible held-out directions.

Eigengene concordance:

- Pearson checks linear agreement of age-aligned module trajectories.
- Spearman checks rank/shape agreement and is less sensitive to amplitude.
- RMSE checks absolute distance between trajectories.
- Strength: good for core-to-core ecological trajectory consistency.
- Weakness: eigengene sign can flip unless explicitly aligned.

Core balance:

- Equal-core downsampling asks whether larger cores dominate module formation.
- Strength: directly targets the R1/R2/ST8/ST13 imbalance concern.
- Weakness: current implementation is one balanced draw, not repeated balanced bootstraps.

Grey burden/module balance:

- Grey fraction and non-grey module count help avoid over-conservative networks where many taxa are dropped into grey.
- Strength: practically important for downstream biology.
- Weakness: lower grey is not automatically better; over-splitting can produce unstable small modules.

## What Is Still Missing

### 1. kME/module membership quality

We should ask whether taxa assigned to a module are actually strongly associated with their module eigengene.

Useful metrics:

- Median assigned kME per module.
- p05 assigned kME per module.
- Fraction of taxa with assigned kME `< 0.2`.
- Fraction of taxa with negative assigned kME.
- Fraction of taxa whose assigned module is also their max-kME module.
- Hub taxa stability: are top-kME taxa recovered under bootstrap?

Why this matters:

- A module can be preserved and have reasonable overlap while still containing many weakly attached taxa.
- kME is the WGCNA-native way to assess fuzzy module membership.
- WGCNA has consensus kME tools, and the manual describes consensus kME as eigengene-based intramodular connectivity/fuzzy membership.

### 2. Eigengene sign alignment

Module eigengenes are first principal components. PC signs are arbitrary.

We already compare eigengene trajectories, but before interpreting direction biologically we should align signs to a reference, probably the training pooled eigengene or GeoB25202_R1.

Useful output:

- Per-module sign used for each core.
- Correlations before and after sign alignment.
- Whether sign flips occur in bootstrap or validation cores.

Why this matters:

- If a module eigengene flips sign, correlation magnitude can still look fine, but biological interpretation of "high" vs "low" module activity becomes inverted.

### 3. Topology quality scoring

We visualize graphs, but the current ranking does not include topology.

Useful metrics:

- Within-module TOM mean/median vs between-module TOM mean/median.
- TOM separation ratio: within / between.
- TOM silhouette-like score by module.
- Modularity of the thresholded TOM graph.
- Fraction of plotted top-TOM edges that are within module.
- Threshold sensitivity at top `0.25%`, `0.5%`, `1%`, `2%`.
- Largest component size and isolate fraction by threshold.

Why this matters:

- A biologically useful module solution should have stronger internal topology than external topology.
- It should not depend on a single plotting threshold.
- The full WGCNA module call uses weighted TOM, but plots and graph summaries require thresholding; threshold sensitivity keeps us honest.

### 4. Preprocessing/depth sensitivity

This is the most important missing safety check.

Compare `exp3`, `exp4`, and `opt5` under a few input transforms:

- Current input: taxon-centered log raw counts with pseudocount.
- Standard sample-wise CLR on raw counts.
- DESeq2/length-normalized log counts.
- Optional: residualize log total reads from each taxon after log transform.

For each transform, record:

- PC1/PC2 correlation with log total reads.
- Module overlap with current `exp3`.
- Grey fraction and module count.
- Preservation count and Z summaries.
- Eigengene concordance.
- kME quality.
- Basic topology separation.

Why this matters:

- If `exp3` is stable across transforms, the current module decision is much stronger.
- If the network changes radically under sample-wise CLR or depth residualization, we should not lock the current modules for downstream biology yet.

### 5. Robust-correlation sensitivity

Current WGCNA uses Pearson. WGCNA guidance often recommends `bicor` when outliers are plausible, with `maxPOutliers = 0.05` or `0.10`.

Given ancient DNA damage and rough count profiles, a limited sensitivity check is reasonable:

- Run `exp3` under Pearson and `bicor`.
- Compare module overlap, grey burden, preservation, and kME.

I would not start with a full bicor parameter sweep unless the focused comparison shows a large difference.

## Reasonable Next Step

The preprocessing/depth sensitivity check has now been moved into the isolated `networkQC/input_evaluation` lane. It showed that depth-related structure remains visible across candidate inputs, so we will keep that as a project-level caveat and move to the next WGCNA-native QC layer.

Next focused comparison:

`baseline` vs `opt5` vs `exp4` vs `exp3`

Purpose:

- `baseline`: original/main-pipeline reference.
- `opt5`: conservative 5-module-style optimized setting.
- `exp4`: runner-up with strongest preservation count.
- `exp3`: current best full-eval setting.

The next QC layer should test whether these settings form strong module memberships and coherent TOM topology, not only stable module labels.

Detailed strategy:

`networkQC/KME_TOPOLOGY_QC_STRATEGY.md`

Add one focused NetworkQC safety script rather than rerunning every sweep:

`networkQC/scripts/10_kme_topology_review.R`

Suggested outputs:

- `networkQC/results/tables/kme_module_membership_summary.tsv`
- `networkQC/results/tables/kme_taxon_membership.tsv`
- `networkQC/results/tables/topology_quality_summary.tsv`
- `networkQC/results/tables/topology_threshold_sensitivity.tsv`
- `networkQC/results/tables/final_qc_integrated_ranking.tsv`
- `networkQC/results/figures/kme_topology_setting_heatmap.png`
- `networkQC/results/KME_TOPOLOGY_QC_REPORT.md`

Scope:

- Run focused checks on `baseline`, `opt5`, `exp4`, and `exp3`.
- Recompute TOM only for the two unique powers needed: 12 and 20.
- Evaluate kME/module membership, grey-rescuable taxa, within/between TOM separation, and threshold sensitivity.
- Avoid a new full parameter sweep unless kME/topology contradict the current `exp3` decision.
- Reuse existing outputs where possible.

## Practical Recommendation

Do not change downstream scripts yet.

Current decision:

- `exp3` is the best parameter set among the evaluated WGCNA configurations.

Safety decision:

- Before locking `exp3` into the main pipeline, run the focused kME/topology review across `baseline`, `opt5`, `exp4`, and `exp3`.

Interpretation stance:

- If `exp3` also holds under kME and TOM topology checks, it becomes the best available consensus network under the current input strategy.
- If `exp4` improves membership/topology enough to offset its lower bootstrap/core-balance scores, consider switching to `exp4`.
- If both expanded settings show weak membership or poor TOM separation, hold `opt5` as the conservative fallback.
- Do not return to `baseline` unless all optimized settings fail, because its grey burden remains too high.

## Notes From the Reasoning Process

The review started from the assumption that our recent QC might already be sufficient. The full-eval results do support that interpretation in terms of the tested parameter space: `exp3` is consistently better than the original baseline and better balanced than the older 5-module options.

The first concern was whether our metrics were too overlap-heavy. That was partly true, but less severe after adding preservation and age-aligned eigengene correlation. The bigger missing WGCNA-native metric is kME/module membership. This matters because WGCNA modules are not only hard clusters; taxa can belong more or less strongly to their assigned module.

The next concern was topology. The graph cleanup showed that the earlier strange-looking plots were largely caused by plotting all isolates at a very sparse TOM threshold. That explains the visual artifact, but it does not replace quantitative topology checks. We still need within/between TOM separation and threshold sensitivity.

The most important concern emerged from reading `01_data_prep.R`: the saved "CLR" matrix is not sample-wise CLR. After checking the current matrix, PC1 showed a strong relationship with total read depth. Since WGCNA correlations are sensitive to sample-wise scaling effects, this is the most important safety issue.

The recommended response is not to discard the QC work. Instead, use the current QC result (`exp3`) as the candidate network and run focused sensitivity checks around it. That gives us a reasonable path between overbuilding and ignoring a real preprocessing risk.
