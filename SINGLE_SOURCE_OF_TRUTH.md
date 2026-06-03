# ROCS Single Source of Truth

Last updated: 2026-05-28 (UTC)
Scope: authoritative project status for `/src` workflow and QC state

## 1) Purpose

This document is the canonical, current-state reference for the ROCS analysis workspace in `/src`.

It replaces conflicting/stale narrative fragments across older notes and summarizes:

- what this project is doing now,
- which workflow is current,
- which parameters/settings are accepted,
- what caveats still apply,
- and what should happen next.

If another markdown file disagrees with this one, treat this file as primary unless newer dated evidence is provided.

## 2) Project Goal

Primary scientific question:

How microbial communities in ancient marine sediments change across glacial/interglacial cycles, and which taxa/modules are most linked to metabolic potential, redox strategy, and climate-associated ecosystem behavior.

Core data context:

- aDNA microbial counts and damage-classified profiles
- climate/proxy metadata (age, d18O/MIS, SST, geochemistry)
- focal cores: ST8, ST13, GeoB25202_R1 (training), GeoB25202_R2 (validation)
- main implemented analysis window is concentrated in the ~0–150 ka interval

## 3) Canonical Repo Layout (Working Copy)

Top-level analysis roots in `/src`:

- `scripts/` main pipeline scripts (`01` to `16`)
- `InputQC/` preprocessing and detection/depth sensitivity diagnostics
- `networkQC/` module/network parameter QC and comparative evaluation
- `results/` main pipeline outputs
- `data/` metadata and medium-size inputs

Important note on duplication:

- There are mirrored markdown/report trees under `/src/MainRepo/rocs`.
- Many are byte-identical copies of `/src/InputQC` and `/src/networkQC` docs.
- For active local work in this environment, prefer `/src/...` paths as canonical.

## 4) Current End-to-End Workflow

Operational script order (from `run_all.sh` and `DATA_SUMMARY.md`):

1. `01_data_prep.R`
2. `02_wgcna.R`
3. `03_hmm_states.R`
4. `04_emp.R`
5. `05_tea_vs_emp.R`
6. `06_figures.R`
7. `06b_bryantfigures.R`
8. `07_taxon_importance.R`
9. `07b_taxon_importance_fuzzy.R`
10. `08_network_statistics.R`
11. `09_driver_integration.R`
12. `10_climate_sensitivity.R`
13. `11_state_networks.R`
14. `12_functional_linkage.R`
15. `13_state_transition_network.R`
16. `14_driver_quadrants.R`
17. `15_state_functional_breakdown.R`
18. `16_final_story_visualization.R`

Start-step normalization supports values like `7`, `07`, `07b`, `010` where defined.

## 5) Current Accepted NetworkQC Decision

Authoritative current module-construction choice:

- Setting ID: `exp3`
- soft power: `12`
- deepSplit: `3`
- mergeCutHeight: `0.25`
- minModuleSize: `20`

Evidence basis:

- `networkQC/results/full_eval/FULL_EVAL_REPORT.md`
- `networkQC/results/KME_TOPOLOGY_QC_REPORT.md`
- `networkQC/results/NETWORK_QC_REPORT.md`
- `networkQC/results/NETWORK_QC_DECISION_REPORT.md`

Key rationale summary:

- best integrated score across full eval + kME + topology
- low grey burden relative to baseline
- strong biological preservation/stability balance
- runner-up `exp4`; conservative fallback `opt5`

## 6) InputQC Decision State

Current operating stance:

- Keep current input as operational anchor for continuity.
- Do not claim the depth/preservation confounding issue is solved.

InputQC core finding:

- leading axes remain strongly associated with technical/detection structure under multiple input variants.
- depth residualization removes depth association but substantially alters module structure.
- ALR variants in quick check were destabilizing for this use case (single-reference especially collapsed module structure).

Evidence basis:

- `InputQC/INPUT_QC_DIAGNOSTIC_REPORT.md`
- `InputQC/INPUT_DECISION_REPORT.md`
- `InputQC/LOW_DETECTION_STRUCTURE_REPORT.md`
- `InputQC/ALR_QUICK_CHECK_REPORT.md`
- `InputQC/INPUT_QC_RESEARCH_PLAN.md`

## 7) Sample-Quality and Imbalance Findings

Rarefaction/depth QC:

- most samples classified as adequate depth
- minority are low-diversity but saturated
- very small subset insufficient for clear call

Core/age imbalance:

- ST8 contributes dense coverage and includes enriched low-detection bands beyond sampling expectation in some age windows
- indicates mixed technical + paleo/preservation structure

ST8 low-taxa review:

- low-taxa ST8 subset is real and concentrated in specific age windows
- lower library concentration is a strong contributor but not sole explanation

Evidence basis:

- `InputQC/rarefaction_depth_qc/RAREFACTION_DEPTH_QC_REPORT.md`
- `InputQC/core_age_imbalance/CORE_AGE_IMBALANCE_REPORT.md`
- `InputQC/st8_low_taxa_review/ST8_LOW_TAXA_REVIEW.md`

## 8) HMM / State Status (Current Snapshot)

Current project summaries indicate a hybrid HMM decision process has been used in recent runs, with state selection balancing BIC ambiguity and held-out behavior.

Important guardrail:

- treat state labels and state-count claims as run-dependent outputs.
- always verify against the latest files in `results/hmm/` and `DATA_SUMMARY.md` before making narrative claims.

Known historical drift exists across docs (older K=5 narratives vs newer K=4 summaries in some reports).

## 9) What Is Considered Stale or Secondary

Potentially stale or mixed-era narrative docs:

- `CODEX.md`
- `GEMINI.md`
- older planning prose in ad hoc notes/slides

These remain useful context, but they are not authoritative for current accepted settings.

Use these as truth anchors first:

1. this document
2. `DATA_SUMMARY.md`
3. `networkQC/results/*` latest dated reports
4. `InputQC/*` latest dated reports
5. current result tables in `results/`, `networkQC/results/tables`, `InputQC/.../results/tables`

## 10) Reproducibility and Runtime Notes

- Workspace runtime here is V1 with V2 tools mounted; do not assume V2 home is active context.
- Squeezr proxy is not active in this session; `rtk` and `sqz` are available.
- Token policy remains: targeted reads, narrow commands, avoid broad output dumps.

Large required inputs may be external in collaborator subtree contexts; verify path resolution in configs before full reruns.

## 11) Decision Rules for Future Changes

Any proposal to replace the current anchor input or `exp3` setting should provide:

1. overlap and stability impact,
2. grey/module balance impact,
3. preservation and age-aligned concordance impact,
4. kME coherence impact,
5. topology separation impact,
6. downstream biological consistency (importance, climate sensitivity, state-network outputs).

No single metric is sufficient.

## 12) Recommended Immediate Next Actions

1. Refresh `DATA_SUMMARY.md` from current outputs after next validated run.
2. Keep `exp3` as default module setting unless a full replacement benchmark is passed.
3. Continue InputQC with explicit sensitivity framing (diagnostic/control vs production candidates).
4. Before manuscript-facing claims, enforce one synchronized snapshot pass:
   - module counts,
   - state counts/labels,
   - top driver tables,
   - climate significance summaries,
   - final figures.

## 13) Quick Command Reference

From `/src`:

```bash
bash run_all.sh
bash run_all.sh --start 07b
bash networkQC/run_network_qc.sh
bash InputQC/run_input_evaluation.sh
bash InputQC/run_depth_and_imbalance_qc.sh
```

## 14) Canonical References (Local)

- `README.md`
- `WORKFLOW_DOCUMENTATION.md`
- `DATA_SUMMARY.md`
- `networkQC/WGCNA_QC_REVIEW_NEXT_STEPS.md`
- `networkQC/results/FULL_EVAL_REPORT.md`
- `networkQC/results/KME_TOPOLOGY_QC_REPORT.md`
- `InputQC/INPUT_QC_RESEARCH_PLAN.md`
- `InputQC/INPUT_QC_DIAGNOSTIC_REPORT.md`
- `InputQC/INPUT_DECISION_REPORT.md`

