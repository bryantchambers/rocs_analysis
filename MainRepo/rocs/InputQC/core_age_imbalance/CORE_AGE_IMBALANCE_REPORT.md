# Core/Age Imbalance Report

- Generated: 2026-05-18 10:32:10 CEST
- Progress log: `InputQC/core_age_imbalance/results/core_age_imbalance.log`

## Main interpretation

This analysis checks whether ST8 low-detection clustering is stronger than expected after accounting for age-bin sample availability. A cluster is more concerning when ST8 has more low-detection samples than expected from its share of samples in the same age bin.

ST8 has `4` age bins with low-detection samples enriched beyond sampling expectation and `1` bins where the pattern is mostly expected from sampling dominance.

## Core sampling density

|core|n|age min|age max|age span kyr|samples per 10 kyr|low-detection pct|sequencing-attention pct|
|---|---:|---:|---:|---:|---:|---:|---:|
|GeoB25202_R1|26|0.69|147.21|146.51|1.77|15.4|0.0|
|GeoB25202_R2|26|0.69|147.21|146.51|1.77|7.7|0.0|
|ST13|48|5.66|148.09|142.43|3.37|8.3|0.0|
|ST8|118|1.52|107.68|106.16|11.12|22.9|2.5|

## Focus age bands

|age band|core|n|median reads|median observed taxa|low-detection pct|sequencing-attention pct|
|---|---|---:|---:|---:|---:|---:|
|40-60 kya|GeoB25202_R1|6|308334|398|33.3|0.0|
|40-60 kya|GeoB25202_R2|6|266844|383|0.0|0.0|
|40-60 kya|ST13|10|446790|501|0.0|0.0|
|40-60 kya|ST8|18|195235|374|22.2|5.6|
|90-110 kya|GeoB25202_R1|4|580310|593|0.0|0.0|
|90-110 kya|GeoB25202_R2|4|730409|600|0.0|0.0|
|90-110 kya|ST13|4|832988|710|0.0|0.0|
|90-110 kya|ST8|16|42592|110|93.8|12.5|

## How to read the fairness labels

- `expected_by_sampling`: a core dominates the age bin, so some clustering from that core is expected.
- `enriched_beyond_sampling`: low-detection samples from that core exceed what its age-bin sample share predicts.
- `mixed_or_background`: low-detection exists but is not strongly enriched beyond sampling.
- `too_sparse`: too few samples to make a fair comparison.

## Outputs

- Sample table: `InputQC/core_age_imbalance/results/tables/core_age_sample_table.tsv`
- 10 kyr imbalance summary: `InputQC/core_age_imbalance/results/tables/age10k_core_imbalance_summary.tsv`
- Focus age bands: `InputQC/core_age_imbalance/results/tables/focus_age_band_core_summary.tsv`
- Core density: `InputQC/core_age_imbalance/results/tables/core_age_density_summary.tsv`
- Winnowing summary: `InputQC/core_age_imbalance/results/tables/sample_winnowing_summary.tsv`
- Figures: `InputQC/core_age_imbalance/results/figures/`
