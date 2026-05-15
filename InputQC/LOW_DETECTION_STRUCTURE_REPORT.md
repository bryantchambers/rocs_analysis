# Low-Detection Structure Report

- Generated: 2026-05-15 10:53:58 CEST
- Scope: qualitative/descriptive investigation of the convergent low-depth/low-detection ordination region
- Progress log: `InputQC/results/low_detection_structure.log`

## Main interpretation

The convergent region in the current ordination is best treated as a low-detection axis, not as a PC-threshold removal target.

Using the primary descriptive group `total_reads < 100000 OR detected_taxa < 250`, the low-detection group contains `37` samples across ages 4.75-133.78 kyr.

Current qualitative classification: **mixed technical and paleoenvironmental/core structure**.

This means filtering may be useful as a sensitivity analysis, but it is not yet justified as a production exclusion rule. The samples may include technical weakness, preservation structure, and real paleoenvironmental signal.

## Primary group composition

Core composition for `total_reads < 100000 OR detected_taxa < 250`:

|core|n|
|---|---:|
|ST8|27|
|ST13|4|
|GeoB25202_R1|4|
|GeoB25202_R2|2|

Site composition:

|site_group|n|
|---|---:|
|ST8|27|
|GeoB25202|6|
|ST13|4|

Glacial class composition:

|glacial_class|n|
|---|---:|
|interglacial_like|22|
|glacial_like|15|

## Filter sensitivity framing

Reasonable candidate filters to test later, not apply immediately:

- `detected_taxa >= 150`: removes the most extreme low-richness tail.
- `detected_taxa >= 250`: stronger richness filter, close to the visual PC1>50 group but based on independent QC.
- `total_reads >= 50000`: conservative depth filter.
- `total_reads >= 100000`: stronger depth filter.
- combined filter `total_reads >= 50000 AND detected_taxa >= 150`: balanced first sensitivity candidate.

A PC1 threshold should remain diagnostic only. It should not be used as the filtering rule.

## Association tests

Association tests are descriptive because samples are time-ordered and cores are not independent. Use them to guide plots and sensitivity analyses, not as final inference.

See `InputQC/results/tables/low_detection_association_tests.tsv`.

## Outputs

- Sample table: `InputQC/results/tables/low_detection_sample_table.tsv`
- Group summaries: `InputQC/results/tables/low_detection_group_summary.tsv`
- Core enrichment: `InputQC/results/tables/low_detection_core_enrichment.tsv`
- Site/glacial summary: `InputQC/results/tables/low_detection_site_glacial_summary.tsv`
- Association tests: `InputQC/results/tables/low_detection_association_tests.tsv`
- Figures: `InputQC/results/figures/*low*`, `*detected*`, `*reads*`, `*shannon*`, `*pc1*`
