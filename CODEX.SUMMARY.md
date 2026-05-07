# CODEX Summary

Snapshot date: 2026-05-04

## What I Read

- Project docs: `README.md`, `CODEX.md`, `GEMINI.md`, `WORKFLOW_DOCUMENTATION.md`, `DATA_SUMMARY.md`, `BryantsNotes.md`.
- Literature/data notes: `lit/fuzzyforest.pdf`, `lit/v91i09.pdf`, `old/Lisiecki2005_copy.txt`, `data/sourcetracker/med-biomes-download.txt`.
- Code: `config.R`, `config_B.R`, `config_ant.R`, `run_all.sh`, `scripts/*.R`, and the old module-composition script in `old/`.
- Current outputs under `results/`, including WGCNA, HMM, EMP/TEA, FuzzyForest, network statistics, driver integration, climate sensitivity, and final figures.

## Repository Purpose

This repository is a reproducible R analysis pipeline for ROCS ancient marine sediment DNA. The scientific goal is to identify microbial ecological states and microbial drivers linked to glacial/interglacial climate cycles, carbon sequestration, metabolic potential, and terminal electron acceptor/redox strategies.

The main biological setup is:

- Data type: ancient environmental DNA from marine cores, with geochemical and climate proxy metadata.
- Main cores: `ST8`, `ST13`, `GeoB25202_R1`, `GeoB25202_R2`.
- Training/analysis cores: `ST8`, `ST13`, `GeoB25202_R1`.
- Validation core: `GeoB25202_R2`.
- `ST5` appears in metadata but is mostly excluded from the main analysis.
- Primary time window in the implemented stage-1 pipeline is <= 150 ka, although the broader project narrative discusses 600 ka.

## Current Pipeline Shape

The intended pipeline has three layers.

1. Foundation:
   - `01_data_prep.R`: filters damaged prokaryote/viral reads, builds a CLR matrix saved as `prokaryotes_vst.rds`, and builds a DESeq2 object for downstream EMP work.
   - `02_wgcna.R`: consensus WGCNA across the three training cores, with GeoB R2 preservation testing.
   - `03_hmm_states.R`: HMM ecological state discovery from residualized module eigengene PCA.
   - `04_emp.R`: EMP and SAP metabolic summaries.
   - `05_tea_vs_emp.R`: TEA/OAP, methanogenesis, methane filter, denitrification, microoxic, and sulfate reduction indices.

2. Driver discovery:
   - `07_taxon_importance.R`: kME, state-membership scores, and baseline random forest.
   - `07b_taxon_importance_fuzzy.R`: FuzzyForest using WGCNA modules as correlated feature blocks.
   - `08_network_statistics.R`: TOM-derived graph metrics, Z-P roles, PageRank, betweenness, closeness, bridging centrality, and nodal efficiency.
   - `09_driver_integration.R`: combines FuzzyForest importance, topology, EMP, and TEA/OAP into driver tiers.

3. Synthesis and figures:
   - `10_climate_sensitivity.R`: GLS/CAR1 climate response for candidate drivers.
   - `11_state_networks.R`: state-specific subnetworks by filtering active taxa from a global adjacency.
   - `12_functional_linkage.R`: joins driver tiers, climate response, bridge taxa, EMP, and TEA classes.
   - `13` through `16`: state-transition and narrative visualizations.

## Literature Takeaways

The FuzzyForest literature strongly fits the problem. Random forest variable-importance scores are biased when predictors are highly correlated, which is exactly the WGCNA/module setting here. FuzzyForest first partitions correlated features into modules, then runs recursive feature elimination within modules, then performs a final random forest on the surviving features. That makes it a reasonable driver-selection layer for `p >> n` microbial feature matrices.

The LR04/Lisiecki text file is a benthic delta18O climate stack used as the climate reference. The `med-biomes-download.txt` file is a large accession/sample table for environmental metagenomes, likely used as source-tracking or external-context data rather than as a direct model input in the main scripts.

## Current Output Snapshot

Current generated outputs show:

- Prokaryote matrix: 214 samples x 1797 taxa.
- WGCNA assignments: 1797 taxa total.
- Current module counts:
  - `grey`: 1196
  - `turquoise`: 237
  - `blue`: 138
  - `brown`: 85
  - `yellow`: 76
  - `green`: 65
- HMM numeric states: 5 numeric states are present.
- HMM label counts currently collapse to 4 labels:
  - `G-A`: 70
  - `IG-B`: 53
  - `IG-C`: 35
  - `IG-E`: 56
- FuzzyForest:
  - 50 selected features.
  - test accuracy: 0.9032258.
  - kappa: 0.8673797.
- Network statistics:
  - 1797 taxa.
  - 50 potential keystones.
  - 4 hidden gems.
  - 436 connectors.
  - 27 hubs.
- Current integrated driver tiers:
  - Tier 1 Super-Driver: 3
  - Tier 2 High Potential: 54
  - Tier 3 Predictive Specialist: 43
  - Peripheral: 1697
- Current Tier 1 taxa:
  - `S__3300027847_13`, turquoise, core heterotrophy.
  - `S__GCA_905182885.1`, turquoise, pelagic heterotroph.
  - `S__GCA_905478185.1`, turquoise, particle heterotroph.
- Climate sensitivity output reports 54 candidate taxa with `FDR < 0.05`.
- State-network output currently has only 4 rows because states 3 and 5 share the same `IG-E` label. `G-A` has the highest inter-module ratio in the current state-network table, but only slightly above `IG-B`.

## Important Drift and Risks

The biggest immediate issue is the HMM label logic. `results/hmm/state_fingerprints.tsv` has five numeric states, but states 3 and 5 are both labeled `IG-E`. Downstream scripts mostly group by `label`, so they collapse two distinct HMM states into one biological category. This affects state-network metrics, bridge-taxon summaries, state functional enrichment, and narrative plots.

The project prose says there are 6 stable non-grey modules including red, but current `module_assignments.tsv` has 5 non-grey modules and no red module. This may be an expected newer run, but the docs and current outputs should be reconciled before making claims about module count.

`DATA_SUMMARY.md` was refreshed on 2026-05-04 from current outputs and should be treated as the current quick-reference memory for output shapes, driver tiers, module-count drift, and state-label warnings.

`run_all.sh` now orchestrates the numbered scripts through step 16, including the driver, state-network, functional-linkage, and final story visualization scripts. It accepts numeric starts like `7`, `07`, and `010`, plus suffix starts like `07b`.

The local shell currently has no `R` or `Rscript` on `PATH`. `run_all.sh` now uses `${RSCRIPT:-Rscript}`, so future execution needs either `Rscript` on `PATH` or an explicit `RSCRIPT=/path/to/Rscript` / container invocation.

`config_B.R` and `config_ant.R` include climate paths and palettes used by `06b_bryantfigures.R`, but current `config.R` does not include `CLIMATE` or `PALETTES`. As written, `06b_bryantfigures.R` is likely not reproducible with the active `config.R`.

Some methods need careful interpretation:

- `11_state_networks.R` creates state subnetworks from a static global adjacency by filtering active taxa, not from state-specific inferred edges. This is useful as a first approximation, but it should not be over-interpreted as a fully state-specific interaction network.
- `08_network_statistics.R` calls `vulnerability` nodal efficiency, but it is not currently the drop in global efficiency after node removal. If the intended claim is "node-removal vulnerability", this needs a true removal-delta calculation.
- EMP uses shifted CLR values as positive weights. That avoids negative contributions, but EMP should be described as a transformed-abundance proxy rather than a direct physical flux.
- Current bridge functional enrichment reports `O2` as the top TEA class for all state bridge sets. The redox-bridge/sulfate-reduction narrative may still be true for subsets, but it needs a targeted check before being presented as the dominant pattern.

## Recommended Approach Going Forward

1. Stabilize reproducibility first.
   - Fix the R runtime path or document the required container invocation.
   - Use the expanded `run_all.sh` to cover scripts 01 through 16 once the runtime is available.
   - Ensure every analysis script writes logs under `logs/`.

2. Fix and rerun HMM labeling before further interpretation.
   - Preserve unique numeric state IDs.
   - Generate unique biological labels for all five states, for example `G-A`, `IG-A`, `IG-B`, `IG-C`, `IG-D`.
   - Rerun state-dependent outputs after this fix: `07`, `07b`, `09`, `10`, `11`, `12`, `13`, `15`, and `16` as needed.

3. Reconcile module counts and docs.
   - Decide whether the accepted current module solution is 5 non-grey modules or whether a missing red module indicates a changed WGCNA run.
   - Update `README.md`, `CODEX.md`, `GEMINI.md`, `WORKFLOW_DOCUMENTATION.md`, and `DATA_SUMMARY.md` after the accepted run is fixed.

4. Strengthen driver selection.
   - Keep FuzzyForest as the primary statistical feature-selection method.
   - Add stability checks across repeated train/test splits or bootstraps.
   - Report class-wise performance, not only overall accuracy, because HMM state classes are imbalanced.
   - Integrate kME/state-fingerprint scores as mechanistic support rather than as a replacement for FuzzyForest.

5. Tighten network interpretation.
   - Decide whether the goal is global topology, state-active topology, or differential state-specific topology.
   - If claiming state-specific networks, consider state-wise correlations/TOMs or differential edges, while watching sample-size limits.
   - If claiming vulnerability, implement node-removal efficiency deltas.

6. Strengthen climate/function synthesis.
   - Use GLS/CAR1 as the default climate-sensitivity model, but explicitly define whether `mis` is being treated as delta18O or a stage proxy.
   - Incorporate SST and XRF/geochemical proxies only where sample coverage and alignment are clear.
   - Separate findings for the <= 150 ka implemented window from claims about the full 600 ka project.

7. Refresh the narrative figures after methodological fixes.
   - The current story is plausible: glacial states appear more coordinated, FuzzyForest finds taxonomic state predictors, and topological integration highlights candidate drivers.
   - The manuscript narrative should wait until HMM labels, module counts, and driver-tier summaries are internally consistent.

## Best Immediate Next Tasks

1. Patch `03_hmm_states.R` state-label assignment so all five numeric states receive unique labels.
2. Run `scripts/inspect_data.R` in the proper R runtime to regenerate `DATA_SUMMARY.md` from code rather than the manual shell refresh.
3. Use the expanded `run_all.sh` as the main orchestrator once the correct `Rscript` path/container is available.
4. Fix `config.R` or `06b_bryantfigures.R` so climate references and palettes are defined in the active config.
5. Rerun downstream driver/network summaries after HMM labels are fixed.
