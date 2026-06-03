# Network Method Comparison

Compares the original WGCNA workflow against the balanced-network workflow.

Primary entrypoint:

```bash
Rscript network_method_comparison/scripts/01_build_comparison_tables.R
```

Main outputs:

- `METHOD_COMPARISON_REPORT.md`
- `RESEARCH_NOTE_WGCNA_IMBALANCE.md`
- `tables/method_level_comparison.tsv`
- `tables/setting_level_comparison.tsv`
- `tables/module_distribution_comparison.tsv`
- `tables/module_stability_comparison.tsv`
- `tables/kme_topology_comparison.tsv`
- `tables/downstream_readiness_checklist.tsv`
- `figures/`

