# Choosing Between Original and Balanced WGCNA Networks

## Current Interpretation

The balanced and original WGCNA approaches answer slightly different questions.

The original network asks: "What modules are strongest in the observed dataset as collected?"

The balanced network asks: "What modules persist when each core and age bin contributes comparable information?"

Because the observed dataset is strongly ST8-heavy, the original network can look more coherent partly because it sees more ST8 information. That is useful if the goal is to model the collected sample set, but risky if the goal is to infer modules that generalize across cores and age structure.

## Current Evidence From This Project

| Evidence axis | Original best: exp3 | Balanced best: top3 | Interpretation |
|---|---:|---:|---|
| Training ST8 fraction | 60.85% | 33.33% | Balanced is much fairer. |
| Bootstrap Jaccard | 0.400 | 0.604 | Balanced is more reproducible under resampling. |
| Grey taxa | 28.66% | 48.53% | Original assigns more taxa to biological modules. |
| Non-grey modules | 8 | 5 | Original is more granular. |
| Strong biological preservation modules | 6/8 | 4/5 | Both are acceptable; balanced has a higher fraction, original has more modules. |
| Median assigned kME | 0.703 | 0.706 | Essentially tied. Balanced is not weak by kME. |
| 5th percentile assigned kME | 0.383 | 0.370 | Similar, with original slightly higher at the lower tail. |
| Assigned module is max-kME module | 83.0% | 87.6% | Balanced assignments are at least as internally coherent. |
| TOM separation ratio | 3.53 | 2.61 | Original modules are more sharply separated topologically. |
| TOM silhouette-like score | 0.717 | 0.617 | Original has cleaner module geometry. |

## Are The Balanced Modules Safe?

Balanced top3 is not an unsafe network. Its modules are defensible, but not perfect.

Reasons it is safe enough to carry forward as a biologically conservative network:

- kME is reasonable: median assigned kME is about 0.706, essentially identical to original exp3.
- Negative assigned kME is absent.
- Only about 0.65% of assigned taxa have kME below 0.2.
- Most taxa are assigned to the module where they have their maximum kME.
- All 5 biological modules have at least moderate preservation, and 4 have strong preservation.
- Bootstrap Jaccard is higher than the original setting.

Reasons it should still be treated cautiously:

- Grey burden is high at 48.5%, so it discards many taxa from biological module interpretation.
- TOM separation is lower than original exp3, meaning the module boundaries are less sharply separated.
- The turquoise module remains the weakest balanced top3 module by bootstrap stability.
- Downstream biology has not yet been rerun on the balanced modules.

## Reasonable Ranges

There is no universal cutoff that proves a WGCNA module is biologically true. The literature supports using multiple checks.

Bootstrap Jaccard:

- Values below about 0.5 are often treated as cluster dissolution.
- Values from 0.6 to 0.75 suggest a real pattern, but imperfect membership certainty.
- Values above about 0.75 are stronger, and above about 0.85 are highly stable.
- By this yardstick, balanced top3 overall mean 0.604 is in the "pattern is present, membership still uncertain" zone. Original exp3 mean 0.400 is weak overall, although individual modules such as green and turquoise are better.

Preservation Zsummary:

- Zsummary below 2: little evidence of preservation.
- Zsummary 2 to 10: moderate preservation.
- Zsummary above 10: strong preservation.
- Both original exp3 and balanced top3 pass this at the module-family level because their biological modules are moderate or strong.

kME:

- kME is the correlation between a taxon's abundance/profile and its module eigengene.
- Values near 1 indicate a strong module member; values near 0 indicate weak membership.
- There is no universal biological threshold, but many WGCNA workflows treat high-kME nodes as hubs, often using project-specific cutoffs such as 0.7, 0.8, or 0.9.
- In this project, balanced top3 and original exp3 have nearly identical median kME, so kME does not argue against balanced.

Grey fraction:

- Grey is the unassigned set, not a biological module.
- Lower grey is useful only if those extra assignments are stable and biologically coherent.
- A low grey fraction can be misleading if it forces weak taxa into modules.

TOM separation:

- TOM separation asks whether taxa within modules share stronger network neighborhoods than taxa across modules.
- Higher TOM separation supports cleaner module geometry.
- Original exp3 is better on this axis, so the original modules are more compact and more sharply separated.

## Recommended Decision Rule

Use balanced top3 as the primary network if the manuscript claim is about cross-core, age-aware biology.

Use original exp3 as the primary network if the manuscript claim is about the strongest module structure in the collected dataset, and explicitly caveat ST8 dominance.

The strongest defensible strategy is to present both:

1. Balanced top3 as the fairness-first primary network.
2. Original exp3 as the higher-resolution sensitivity network.
3. Promote only biological conclusions that reproduce across both or are clearly explained as method-specific.
4. Use downstream biology as the final arbiter: HMM state separation, ecological drivers, functional coherence, and whether module conclusions are ST8-dependent.

## Literature Support

- WGCNA FAQ: heterogeneous samples can strongly affect unsupervised WGCNA; categorical sources with enough samples can be handled by consensus/module-per-group analysis, and unwanted large variation may need adjustment.
- WGCNA 1.74 manual: WGCNA implements module detection, module eigengenes, module membership, consensus modules, TOM, and preservation statistics.
- Langfelder and Horvath 2008: WGCNA modules are summarized by eigengenes and module membership/kME; modules are not only clusters, they are eigengene-linked network summaries.
- Langfelder et al. 2011/module preservation literature: preservation should be evaluated with density/connectivity statistics; Zsummary >10 is strong and 2-10 is moderate evidence.
- Cluster bootstrap literature: Jaccard stability is useful, but stability alone does not prove validity; very stable clusters can still be biologically wrong if the input structure is biased.
- Microbial network literature: compositionality, sparsity, environmental covariates, and heterogeneous samples make microbial association networks hard to interpret as direct ecological interactions. WGCNA modules should therefore be interpreted as co-abundance modules unless validated downstream.

## Sources

- WGCNA FAQ: https://edo98811.github.io/WGCNA_official_documentation/faq.html
- WGCNA 1.74 manual: https://cran.r-universe.dev/WGCNA/doc/manual.html
- WGCNA package paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/
- Module validation review: https://www.nature.com/articles/srep15258
- fpc clusterboot guidance: https://rdrr.io/cran/fpc/man/clusterboot.html
- SPIEC-EASI paper: https://arxiv.org/abs/1408.4158
- FlashWeave paper: https://www.sciencedirect.com/science/article/pii/S2405471219302716
- Microbial co-occurrence caution: https://link.springer.com/article/10.1186/s12859-019-2915-1
