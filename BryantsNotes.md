# Bryant's Analysis Notes: Ancient Ocean Ecosystems

## HMM Hybrid Decision (2026-05-05)

We switched `03_hmm_states.R` to a hybrid strategy to favor robust, low-noise state calls in damaged ancient DNA:

1. Fit candidate HMMs (K=2..5) on training cores only (`ST8`, `ST13`, `GeoB25202_R1`).
2. Score held-out `GeoB25202_R2` using fixed trained parameters.
3. Select K from BIC-ambiguous models (ΔBIC <= 10) by held-out log-likelihood per sample and transition stability.
4. Refit selected K on all cores to produce final state labels for downstream analysis.

Current run summary:
- Train-BIC best model: **K=5** (BIC 1185.36)
- Hybrid-selected model: **K=4**
- Held-out logLik/sample: **K=4 = -15.392**, **K=5 = -15.477**
- Held-out confidence (mean max posterior): **0.990** (K=4)
- Held-out switch rate: **4.17 transitions per 100** (K=4)

Interpretation: K=4 gives slightly better holdout predictive behavior with clean transitions and high confidence, while staying within the BIC-ambiguous band of K=5.

## Topological Metrics & Driver Paradoxes

### 1. The "Peripheral Super-Driver" Paradox
**Observation:** Some taxa are identified as "Super-Drivers" (high predictive importance) but are labeled as "Peripheral" in the Z-P (Within-module Degree vs. Participation) plots.

**Explanation:** 
*   **Local vs. Global Influence:** Z-score and P-score are *local* measures. A "Peripheral" node has a low degree relative to its own module. 
*   **Global Influence (PageRank):** PageRank and Betweenness are *global* measures. A taxon can be "Peripheral" locally (not a hub for its module) but still be globally influential if it acts as a critical bridge between modules or is connected to other influential nodes.
*   **Hidden Gems:** These are the "quiet influencers"—taxa that don't dominate their local module but are vital for whole-system stability.

### 2. Tiered Candidate Identification
To avoid being too conservative (only 9 Super-Drivers), we use a tiered system:

*   **Tier 1 (Super-Drivers):** The "Gold Standard." Top 10% in Statistical Importance (Fuzzy Forest) AND Top 10% in Composite Topological Influence (PageRank + Vulnerability + Betweenness).
*   **Tier 2 (High Potential):** Top 20% in both, OR taxa identified as "Keystones" or "Hidden Gems" in the network analysis. These are likely important but might have slightly lower statistical signal or local degree.
*   **Tier 3 (Predictive Specialists):** High Statistical Importance (>90th percentile) but lower topological influence. These are great at predicting HMM states but aren't central to the network structure.

## Integration Logic (Script 09 Update)
*   **Composite Topology:** Instead of just PageRank, we now use the average percentile of PageRank, Vulnerability (Nodal Efficiency), and Betweenness.
*   **Inclusion:** We explicitly pull in the `is_hidden_gem` and `is_potential_keystone` flags from the network statistics.

## Path A: Climate Sensitivity (GLS Analysis)
*   **Metric:** Generalized Least Squares (GLS) with CAR1 autocorrelation.
*   **Discovery:** 54 of our top 100 candidates are statistically sensitive to climate cycles (FDR < 0.05).
*   **Pattern:** Higher metabolic capacity (EMP) often correlates with positive d18O coefficients, suggesting "metabolic resilience" is a key trait for surviving Glacial (cold) periods.

## Path B: State-Specific Sub-networks
*   **Metric:** Inter-module connectivity ratio and bridge taxa identification.
*   **Discovery:** The Glacial state (**G-A**) has the highest inter-module connectivity (Ratio = 0.30).
*   **Bridge Roles:** The **Yellow module** acts as the primary "structural glue" during Glacial periods, bridging the Turquoise "Core" module to other functional units.

## Functional Linkage & Biological Insights
*   **Functional Super-Hubs:** 116 taxa identified as having both high network influence (Integrated Score) and high metabolic versatility (EMP).
*   **The "Redox Bridge":** Glacial states are enriched for bridge taxa carrying anaerobic pathways (Sulfate Reduction/Methanogenesis). 
*   **Climate-Function Coupling:** Microbes with the highest metabolic capacity are disproportionately favored during Glacial cycles, acting as the stable "engines" of the ecosystem when environmental forcing is high.

## kME / Module Membership QC (2026-05-13)

### What kME means

`kME` means module eigengene membership.

A module eigengene is the main shared abundance pattern of a module. In WGCNA terms, it is the first principal component of all taxa assigned to that module.

For each taxon, kME asks a simple question:

"Does this taxon move with the overall pattern of this module?"

Technically:

```text
kME(taxon, module) = cor(taxon abundance profile, module eigengene)
```

So a high kME means the taxon tracks the module closely. A low kME means the taxon was assigned to the module but does not strongly follow it. A negative kME means the taxon moves opposite the module eigengene, which is usually suspicious for an assigned module member.

### What we tested

We compared four WGCNA parameter sets in `networkQC`:

*   **baseline:** original/reference WGCNA setup.
*   **opt5:** conservative optimized 5-module-style setup.
*   **exp4:** runner-up expanded setting.
*   **exp3:** current best expanded setting.

For each setting, we asked:

*   Are biological module members strongly attached to their assigned module?
*   Are many taxa weakly attached?
*   Are any taxa negatively attached?
*   Is the assigned module usually the taxon's strongest kME module?
*   Are many grey taxa actually close to biological modules and possibly "rescuable"?

This is important because overlap, preservation, and bootstrap stability tell us whether modules are reproducible, but kME tells us whether the modules are internally coherent.

### Main kME metrics

*   **bio_median_assigned_kME:** the typical strength of biological module membership. Higher is better.
*   **bio_p05_assigned_kME:** the lower-tail membership strength. This asks whether the weakest 5% of assigned biological members still have reasonable membership.
*   **bio_frac_assigned_kME_lt_0_2:** the fraction of biological taxa with weak assigned membership. Lower is better.
*   **bio_frac_negative_assigned_kME:** the fraction of assigned biological taxa moving opposite their module. Lower is better.
*   **bio_frac_assigned_is_max_kME:** the fraction of taxa whose assigned module is also their strongest kME module. Higher is better.
*   **bio_median_kME_margin:** assigned kME minus next-best biological kME. Higher means the taxon is less ambiguous between modules.
*   **grey_rescuable_fraction:** the fraction of grey taxa with decent kME to a biological module. High values suggest the setting may be throwing biologically useful taxa into grey.

### Outcomes by parameter set

**baseline**

The biological modules look strong by kME, but this is misleading in context. Biological median kME was about `0.783`, and assigned-is-max was about `0.897`. However, baseline put `66.56%` of taxa into grey, and `37.9%` of grey taxa looked potentially rescuable. Interpretation: baseline gets clean biological modules partly by excluding a large amount of usable signal into grey.

**opt5**

This is conservative and acceptable. Biological median kME was `0.689`, assigned-is-max was about `0.838`, and there was no weak or negative biological kME problem. Grey burden was lower than baseline but still `34.06%`. Interpretation: opt5 remains a reasonable fallback if expanded settings look over-fragmented, but the expanded settings did not fail.

**exp4**

This has good kME but is less decisive than exp3 overall. Biological median kME was `0.710`, slightly higher than exp3, but assigned-is-max was lower at `0.792`, and kME margin was smaller. Interpretation: members are attached, but some taxa are more ambiguous between modules. Exp4 has strong preservation, but its topology score was weaker than exp3.

**exp3**

This has the best overall balance. Biological median kME was `0.703`, p05 assigned kME was `0.383`, there were no weak or negative biological kME issues, assigned-is-max was `0.830`, and grey-rescuable fraction was only `0.181`. Interpretation: exp3 does not have the highest single kME number, but it has strong membership while also keeping grey burden low and topology strong.

### Why this supports exp3

The kME results show that exp3 is not simply creating extra modules by over-splitting noise. Its biological modules have coherent members, few ambiguous assignments, no weak or negative membership flags, and fewer potentially useful taxa trapped in grey than baseline.

This matters for later steps because:

*   Fuzzy Forest module-aware feature selection becomes more trustworthy.
*   Module eigengene fingerprints for HMM states are more interpretable.
*   kME-weighted driver scoring is meaningful because assigned taxa actually track their module.
*   Hub taxa and bridge taxa are less likely to be artifacts of loose module assignment.
*   Functional enrichment by module has a stronger biological basis.

Current practical conclusion: `exp3` remains the best current module set. The sequencing-depth artifact remains a standing caveat, but kME does not argue against exp3. It strengthens the case for using it.

---
*Last Updated: 2026-05-13*
