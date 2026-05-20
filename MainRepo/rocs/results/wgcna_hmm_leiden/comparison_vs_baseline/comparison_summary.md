# Baseline vs Leiden alternative workflow comparison

Generated: 2026-05-06 21:47:55

Leiden modules here are a graph-partitioning sensitivity analysis and are not directly equivalent
to dynamic-tree modules from consensus WGCNA.

## Grey burden
- Baseline grey: 1196 / 1797 (66.56%)
- Leiden grey: 148 / 1797 (8.24%)
- Delta grey fraction (Leiden - baseline): -58.32 percentage points

## Key downstream metrics
- hmm_states_main: baseline=5; leiden=5; delta=0
- loco_stable_n: baseline=10; leiden=11; delta=1
- n_modules_non_grey: baseline=5; leiden=6; delta=1
- n_taxa_main: baseline=1797; leiden=1797; delta=0
- r1_r2_state_pct_match: baseline=92.3076923076923; leiden=84.6153846153846; delta=-7.69230769230771
- tea_module_enrichment_calls_n: baseline=24; leiden=23; delta=-1
- xrf_significant_n: baseline=18; leiden=24; delta=6

## Interpretation
- Check whether grey burden is reduced while preserving high-level downstream signals.
- Strong shifts in HMM/XRF/TEA metrics indicate sensitivity to module-construction choice.
