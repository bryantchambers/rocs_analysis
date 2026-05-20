# WGCNA+HMM clean workflow execution report

Generated: 2026-05-06 21:47:54

## Data scope
- Main-window samples (<= 150 kyr): 214
- Taxa retained after prevalence filter: 1797
- Non-grey modules: 6

## Module preservation (GeoB25202_R1 -> GeoB25202_R2)
- Strong: 7
- Moderate: 1
- Weak: 0
- Median Zsummary: 14.787

## HMM states
- Operational harmonization mode for states/downstream analyses: corewise_z
- Number of states in main analysis: 5
- R1/R2 interpolated concordance (% match): 84.62
- Main-window state fingerprints file present: yes (state_fingerprints_main.tsv)

## XRF and LOCO
- Significant module-XRF associations (q < 0.10): 24
- LOCO-stable associations: 11

## Extended-time projection
- Samples in extended projection: 415
- Samples older than 150 kyr: 201
- Locked-overlap policy: states at <= 150 kyr are inherited from main analysis; only > 150 kyr states are newly inferred
- Inherited main-window states (state_origin=inherited_main_window): 214
- Newly inferred older states (state_origin=inferred_extended_only): 201
- Overlap final-state agreement with main (after lock): 100.00%
- Selected harmonization mode for extended inference: corewise_z
- Cores requiring scaling fallback (no core-specific main-window rows): 0
- Mode evaluation table: extended_harmonization_mode_evaluation.tsv
- Locked-overlap diagnostics (main vs inferred vs final): extended_locked_overlap_diagnostics.tsv
- State-origin summary (inherited vs inferred): extended_state_origin_summary.tsv
- Scaling parameters used in projection: extended_corewise_z_scaling_parameters_used.tsv
- Scaling fallback diagnostics by core: extended_corewise_z_scaling_fallback_diagnostics.tsv
- Older-period diagnostics by core: extended_older_period_diagnostics.tsv

## Output index
- Main results: results/wgcna_hmm/main/
- Main HMM state fingerprints: results/wgcna_hmm/main/state_fingerprints_main.tsv
- Main TEA results: results/wgcna_hmm/main/tea/
- Taxon-level TEA traits: results/wgcna_hmm/main/tea/tea_taxon_metabolic_calls.tsv
- Extended results: results/wgcna_hmm/extended/
- Report tables: results/wgcna_hmm/reports/
