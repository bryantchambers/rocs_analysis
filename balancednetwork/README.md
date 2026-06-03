# balancednetwork

Isolated experimental workflow for site-age balanced WGCNA construction.

## Goal

Build and evaluate a balanced consensus network so ST8 does not dominate
training sample composition by count or by age-bin representation.

## Design

- Training cores: `ST8`, `ST13`, `GeoB25202_R1`
- Validation core: `GeoB25202_R2` (held out)
- Age-bin balancing: default `10` kyr bins
- Quota per bin: minimum sample count across training cores in that bin
- Baseline balanced build: sample without replacement
- Bootstrap balanced stability: sample with replacement inside each core-bin

## WGCNA settings

Fixed to `exp3`:

- `soft_power = 12`
- `deepSplit = 3`
- `mergeCutHeight = 0.25`
- `minModuleSize = 20`

## Run

From `/src`:

```bash
bash balancednetwork/run_balancednetwork.sh --mode=build
bash balancednetwork/run_balancednetwork.sh --mode=final
bash balancednetwork/run_balancednetwork_qc.sh
bash balancednetwork/run_balancednetwork_full_eval.sh --mode=build
bash balancednetwork/run_balancednetwork_full_eval.sh --mode=final
```

Optional overrides:

- `N_BOOT` (default build: `100`, final: `1000`)
- `N_PERM` (default build: `200`, final: `700`)
- `BIN_WIDTH_KYR` (default: `10`)
- `FORCE=1` (rerun and overwrite existing outputs)
- `TOP_N` (full eval candidate count; default `5`)
- `N_BOOT` default in `run_balancednetwork_full_eval.sh`: build `100`, final `500` (initial run target)

## Outputs

- `balancednetwork/results/tables/` sample balancing design tables
- `balancednetwork/results/wgcna/` balanced consensus outputs
- `balancednetwork/results/stability/` bootstrap stability outputs
- `balancednetwork/results/reports/BALANCED_NETWORK_QC_REPORT.md`
- `balancednetwork/results/qc/tables/` sweep + decision matrix tables
- `balancednetwork/results/qc/BALANCED_QC_DECISION_REPORT.md`
- `balancednetwork/results/qc/full_eval/` per-candidate full evaluation outputs
- `balancednetwork/results/qc/BALANCED_FULL_EVAL_REPORT.md`
