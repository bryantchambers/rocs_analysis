# WGCNA Imbalance Research Note

## Question

Can unequal sample density, especially ST8 dominance, make the original WGCNA network look more stable than a more sampling-balanced network?

## Method Background

- WGCNA builds weighted correlation networks and summarizes modules with eigengenes, kME/module membership, and topological overlap. Core references: Langfelder and Horvath 2008 (<https://link.springer.com/article/10.1186/1471-2105-9-559>) and the WGCNA manual (<https://cran.r-universe.dev/WGCNA/doc/manual.html>).
- The WGCNA FAQ warns that strong categorical drivers, batch effects, or biological heterogeneity can dominate correlations; it recommends inspecting sample clustering and considering adjustment or consensus-style analyses when heterogeneity is strong (<https://edo98811.github.io/WGCNA_official_documentation/faq.html>).
- Consensus WGCNA is designed to find modules that recur across datasets or groups. In this project, cores/age bins behave like structured groups, so a balance-aware comparison is more defensible than treating raw sample count as neutral information.

## Project-Specific Risk

- Original training sample counts are imbalanced: ST8 fraction is 60.8%; balanced training ST8 fraction is 33.3%.
- InputQC already flags ST8 as much more densely sampled and shows a 90-110 kya band where ST8 low-detection samples are enriched beyond sampling expectation.
- Therefore, original-network stability can mean either true ecological robustness or repeated recovery of an ST8-heavy correlation structure.

## Practical Interpretation Rule

Treat bootstrap stability as necessary but not sufficient. A fair network should retain stability after core/age balancing, show reasonable grey burden, preserve modules in GeoB25202_R2, maintain kME/TOM coherence, and support downstream ecological state/driver interpretations.
