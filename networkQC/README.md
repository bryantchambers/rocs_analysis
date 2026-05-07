# networkQC

Isolated workflow for network/module quality control and method comparison.

## Scope

This workflow compares two module construction strategies built on the same expression inputs:

1. Consensus WGCNA module assignment.
2. Leiden clustering on a WGCNA-style adjacency/TOM-derived network.

It also runs foundational QC:

- soft-threshold behavior,
- sample and gene clustering diagnostics,
- dynamic tree cut sensitivity,
- per-core eigengene agreement,
- preservation summaries,
- parameter sweep summaries.

## Run

From `/src`:

```bash
bash networkQC/run_network_qc.sh
```

Outputs are written to:

- `networkQC/results/tables`
- `networkQC/results/figures`

