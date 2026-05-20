# Broad deep-sediment exclusion sensitivity comparison

## Exclusion definition used
- curated diagenetic rule:
  - `display_category == "Diagenetic"`
  - `signal_source %in% c("diagenetic_in_situ", "benthic_resident")`
  - `confidence_score >= 0.85`
- plus extra excluded phyla: `p__Aerophobota`, `p__Atribacterota`, `p__Asgardarchaeota`

## Exclusion counts
- Excluded taxa (unique): 472
- Excluded manifest rows (taxon x reason): 472
- Excluded taxa that would pass prevalence in main window: 200

## Key metric deltas (sensitivity - baseline)
- hmm_states_main: baseline=5, sensitivity=5, delta=0
- loco_stable_n: baseline=10, sensitivity=8, delta=-2
- n_modules_non_grey: baseline=5, sensitivity=4, delta=-1
- n_taxa_main: baseline=1797, sensitivity=1559, delta=-238
- r1_r2_state_pct_match: baseline=92.3076923076923, sensitivity=84.6153846153846, delta=-7.69230769230771
- tea_module_enrichment_calls_n: baseline=24, sensitivity=19, delta=-5
- xrf_significant_n: baseline=18, sensitivity=14, delta=-4

## Module membership stability
- Shared taxa compared: 1559
- Taxa changing assigned module: 193 (12.38%)

## HMM state comparison
- Exact state label match across samples: 0.00%
- Label-invariant pairwise partition agreement: 86.25%

## Brown-associated TEA pattern check
- brown OAP include_grey: baseline=depleted (q=3.432645e-19), sensitivity=enriched (q=1.493314e-05)

## Interpretation (conservative)
- Architecture-level metrics are compared descriptively; directional TEA changes indicate sensitivity for specific module-index associations.
- These diagnostics are not causal attribution and should be interpreted with classification uncertainty in mind.
