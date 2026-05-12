# QC Notes Cheat Sheet (WGCNA + NetworkQC)

Last updated: 2026-05-08

## Why this exists

This is a quick-reference note file for interpreting our `networkQC` outputs and choosing robust module settings before downstream biology.

---

## 1) What the Network QC report is telling us (this run)

Source: `networkQC/results/NETWORK_QC_REPORT.md`

- Training samples: `189`
- Taxa: `1797`
- Non-grey WGCNA modules: `5`
- Grey fraction: `66.6%` (high)
- Sweep tested: `72` parameter combinations
- Best (lowest grey) sweep candidate reached grey ~`28.7%`

Interpretation:
- The current baseline WGCNA settings produce many unassigned taxa (large grey module).
- Parameter tuning can substantially reduce grey burden.

---

## 2) Core terms to remember

### KME / Module membership
Correlation of each taxon with a module eigengene. High absolute kME means the taxon is strongly associated with that module profile.

### Module eigengene (ME)
The first principal component of taxa in a module. Think of it as the module’s summary time-series.

### Grey module
In WGCNA, `grey` means “not confidently assigned to a non-grey module.”  
Important: grey can still contain biological signal; it is not automatically “noise only.”

### Consensus WGCNA
Finds modules that are reproducible across multiple datasets/cores. More conservative than single-dataset clustering.

---

## 3) Metrics you asked about (plain language + formula intuition)

### Jaccard overlap
How similar two module member sets are:

`Jaccard = |A ∩ B| / |A ∪ B|`

- `1.0` = identical membership
- `0.0` = no overlap

Used in bootstrap stability to ask: “Do we recover the same module repeatedly?”

### p05 and p95 (in stability table)
- `p05`: 5th percentile of Jaccard across bootstrap reruns
- `p95`: 95th percentile

Meaning:
- low `p05` near 0 => module frequently unstable/collapsing in some reruns
- narrow range with high median => stable module definition

### baseline_size
Number of taxa in the module in the baseline (reference) run.

### bootstrap_size_median
Median size of the best-matched module across bootstrap reruns.

If very different from baseline_size, module boundaries may be unstable.

---

## 4) Concordance metrics (age-aligned eigengenes)

### Pearson r
Linear agreement of aligned ME trajectories.
- High positive => move together linearly

### Spearman rho
Rank/shape agreement (less sensitive to scaling and some nonlinear behavior).
- High positive => same monotonic trend

### RMSE
Average magnitude gap between aligned trajectories.
- Lower => closer absolute values

Why all three:
- Pearson: linear co-movement
- Spearman: trend/rank similarity
- RMSE: absolute mismatch size

---

## 5) How age alignment works (in our scripts)

For each module ME and pair of cores:
1. Convert sample depths to `age_kyr`.
2. Find overlapping age window.
3. Build a common age grid (100 points).
4. Interpolate each core ME onto that same grid.
5. Compute Pearson/Spearman/RMSE on paired grid points.

This avoids unfair direct comparison of unevenly sampled cores.

---

## 6) Interpreting this run’s concordance

From `NETWORK_QC_REPORT.md`:

- Strongest agreement: `GeoB25202_R1` vs `GeoB25202_R2` (Pearson ~0.85)
- Moderate: `ST13` vs GeoB cores
- Weak: `ST8` vs others

Interpretation:
- ST8 may represent a distinct structure or noisier profile relative to others.
- Consensus constraints could push more taxa into grey if ST8 differs strongly.

---

## 7) Soft-threshold and dynamic tree cut refresher

### Soft-threshold power
Controls how strongly correlations are emphasized in adjacency:

`adj ~ |cor|^power` (signed variant in our setup)

High power:
- suppresses weak correlations
- can increase sparsity and fragmentation
- can increase grey burden if over-strong

Low power:
- keeps more weak edges
- may blur module boundaries

### Dynamic tree cut controls
- `deepSplit`: higher => more/smaller modules (more splitting)
- `mergeCutHeight`: higher => more merging after split
- `minModuleSize`: minimum allowed module size

Need to tune these jointly with power.

---

## 8) Why our decision matrix preferred power=12

From `NETWORK_QC_DECISION_REPORT.md`, top candidates cluster at:
- `power=12`, `deepSplit=1`, `mergeCutHeight=0.10–0.15`, `minModuleSize=20–30`

Reason:
- much lower grey burden than baseline
- module count near target (5 non-grey)
- reasonable module size behavior

This is a candidate set, not final truth. Must pass stability + downstream checks.

---

## 9) Known caveats from current `networkQC` implementation

1. Parameter sweep summary does **not yet** include full bootstrap stability per sweep setting.
2. Decision matrix currently uses global concordance benchmark, not per-setting concordance.
3. Leiden implementation currently yields many tiny communities in this setup; needs stronger filtering/tuning for fair biological comparison.

Action:
- run full stability and downstream checks on top 2–3 WGCNA settings before locking one.

---

## 10) Quick “what to do next” checklist

1. Re-run WGCNA + 02b stability on top 3 parameter sets.
2. Compare:
   - median Jaccard (target higher),
   - p05 floor (avoid frequent collapse),
   - biological preservation classes,
   - age-aligned concordance.
3. Run 07b/09/10/11/12 for each finalist.
4. Choose setting with best joint technical + biological robustness.

---

## 11) Questions/confusions explicitly tracked for follow-up

These were raised and should stay visible:

- “What do Jaccard / p05 / p95 mean?”
- “Why are baseline_size and bootstrap_size_median different?”
- “Does large grey mean no biology?”
- “If grey has predictive power in FuzzyForest, how should we interpret it?”
- “How do Pearson/Spearman/RMSE differ and why use all?”
- “How does age-grid alignment actually work?”
- “Could consensus pressure across cores push taxa into grey?”
- “Why did the full graph plots look like a random shotgun pattern?”
- “Could the original top-5 full evaluation miss better settings with more than 5 modules?”

---

## 12) Network graph plotting lesson

The original `full_graph_*.png` plots looked like random scatter because the graph included all `1797` taxa, but at the plotted top-edge threshold most nodes were isolates.

Current diagnostics at the top 0.5% TOM edge threshold:

- Baseline: 455 non-isolated nodes, 1342 isolates
- opt5: 349 non-isolated nodes, 1448 isolates
- Leiden: 349 non-isolated nodes, 1448 isolates

Interpretation:
- The layout algorithm spreads isolated nodes around the canvas.
- Those isolated nodes visually swamp the real connected components.
- The network is not necessarily wrong; the all-node visualization was misleading.

Better plots:
- remove isolates,
- show largest component separately,
- show non-grey non-isolates separately,
- show module-level meta-graph.

Generated corrected plots:
- `networkQC/results/figures/layout_diagnostics/*_nonisolates_fr.png`
- `networkQC/results/figures/layout_diagnostics/*_largest_component_fr.png`
- `networkQC/results/figures/layout_diagnostics/*_module_meta_graph.png`

---

## 13) Top-N parameter sweep caveat

The first full evaluation used the decision-matrix top 5. Those were all in the same parameter region:

- `power=12`
- `deepSplit=1`
- `5` non-grey modules
- grey fraction about `34%`

This was reasonable if we strongly prefer a 5-module answer. But it may miss better solutions if we allow 6-10 modules.

Promising excluded settings:

- `power=12, deepSplit=2, mergeCutHeight=0.25, minModuleSize=30`: 6 modules, grey about `33.2%`
- `power=12, deepSplit=2, mergeCutHeight=0.25, minModuleSize=20`: 7 modules, grey about `29.1%`
- `power=12, deepSplit=3, mergeCutHeight=0.25, minModuleSize=20`: 8 modules, grey about `28.7%`
- `power=12, deepSplit=3, mergeCutHeight=0.20, minModuleSize=20`: 9 modules, grey about `28.7%`

Next best check:
- run a second full evaluation on a diverse candidate shortlist, not just the top 5 under the 5-module prior.

---

## 14) TOM edge threshold: plotting versus module construction

The `top 0.5% TOM edge threshold` is a graph-plotting choice, not a WGCNA
module-construction parameter.

Where it appears:
- `networkQC/scripts/08_graph_layout_diagnostics.R` computes adjacency, then TOM,
  then keeps only edges above the 99.5th percentile of TOM values for plotting.

Where module construction happens:
- `scripts/02_wgcna.R` and `networkQC/scripts/02_wgcna_parameter_sweep.R` use
  `blockwiseConsensusModules(...)`.
- That function uses the full weighted network structure to build modules. It is
  not limited to the top 0.5% plotted TOM edges.

Important distinction:
- `power` changes the adjacency matrix by raising correlations to a soft-threshold
  power.
- WGCNA converts adjacency into TOM, which measures shared-neighborhood similarity.
- WGCNA clusters taxa using TOM-based dissimilarity, usually `1 - TOM`.
- `deepSplit` controls how aggressively the dendrogram is cut into modules.
- `minModuleSize` prevents tiny modules.
- `mergeCutHeight` merges modules with similar eigengenes after the initial tree cut.

So:
- The optimized WGCNA parameters affect module identities.
- The top-edge TOM threshold affects how much of the network we can see clearly in
  a static plot.

Cheat-sheet definitions:
- Correlation matrix: pairwise similarity between two taxa across samples/time.
- Adjacency matrix: softened network strength after applying the WGCNA power. In
  a signed network, strong positive correlations become strong edges; weak or
  negative relationships become weak edges.
- TOM, or topological overlap matrix: neighborhood-aware similarity. Two taxa have
  high TOM if they are directly connected and/or connected to many of the same
  other taxa.

Why plots need filtering:
- With about `1797` taxa, a full weighted graph has more than 1.6 million possible
  pairwise edges.
- Plotting all edges is unreadable.
- Keeping the strongest TOM edges makes the graph visible, but it is a diagnostic
  simplification of the fitted network, not the network used to call modules.

How to interpret graph plots:
- Use them to ask whether strong edges mostly stay within modules.
- Use them to spot isolates, bridges, and strange module mixing.
- Do not judge module validity from one edge threshold alone.
- Compare several TOM thresholds, non-isolate plots, largest-component plots, and
  module-level meta-graphs.
