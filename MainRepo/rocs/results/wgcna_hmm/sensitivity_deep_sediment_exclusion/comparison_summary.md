# Deep-sediment exclusion sensitivity comparison

## Exclusion definition used
- `display_category == "Diagenetic"`
- `signal_source %in% c("diagenetic_in_situ", "benthic_resident")`
- `confidence_score >= 0.85`

## Exclusion counts
- Candidate taxa meeting criteria in classification table: 425
- Excluded taxa present in selected samples: 425
- Excluded taxa that would otherwise pass prevalence in main window: 200

## Key metric deltas (sensitivity - baseline)
- n_taxa_main: baseline=1797.0, sensitivity=1597.0, delta=-200.0
- n_modules_non_grey: baseline=5.0, sensitivity=5.0, delta=0.0
- hmm_states_main: baseline=5.0, sensitivity=5.0, delta=0.0
- r1_r2_state_pct_match: baseline=92.3076923076923, sensitivity=96.1538461538462, delta=3.8461538461538964
- xrf_significant_n: baseline=18.0, sensitivity=18.0, delta=0.0
- loco_stable_n: baseline=10.0, sensitivity=11.0, delta=1.0
- tea_module_enrichment_calls_n: baseline=24.0, sensitivity=22.0, delta=-2.0

## Module membership stability
- Shared taxa compared: 1597
- Taxa changing assigned module: 174 (10.90%)

## HMM state comparison
- Exact state label match across samples: 0.47%
- Label-invariant pairwise partition agreement: 97.53%
- Note: direct state label IDs may be permuted between runs.

## Brown-associated TEA pattern check
- brown OAP include_grey: baseline=depleted (q=3.43264505117695e-19), sensitivity=enriched (q=6.297432879717059e-16).

## Interpretation (conservative)
- Core architecture (non-grey module count, state count) is similar, but at least one key brown-associated TEA signal changes direction, indicating material sensitivity for that specific pattern.
- These results are descriptive sensitivity diagnostics and do not by themselves identify causality.