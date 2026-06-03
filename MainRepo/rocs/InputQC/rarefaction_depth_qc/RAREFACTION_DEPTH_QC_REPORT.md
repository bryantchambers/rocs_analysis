# Rarefaction and Sequencing Depth QC Report

- Generated: 2026-05-18 10:32:02 CEST
- Progress log: `InputQC/rarefaction_depth_qc/results/rarefaction_depth_qc.log`
- Packages: `phyloseq`, `vegan`, `iNEXT`

## Main interpretation

This analysis asks whether low-read / low-taxa samples appear undersequenced, or whether they look like genuinely low-diversity samples that have mostly saturated at their current depth.

Important caveat: Chao1 is weak for this processed count matrix because many samples have no singleton taxa; when Chao1 equals observed richness, the report relies more on rarefaction slope and observed depth/richness behavior.

## QC class summary

|QC class|n|pct|median reads|median observed taxa|median slope per 10k|
|---|---:|---:|---:|---:|---:|
|adequate_depth|177|82.7|511586|541|0.00|
|low_diversity_but_saturated|34|15.9|41495|108|0.00|
|insufficient_for_call|3|1.4|6403|20|0.00|

## Threshold sensitivity summary

These thresholds should be treated as sensitivity candidates, not automatic production exclusions.

|threshold|flag value|n|pct|median reads|median observed taxa|pct likely undersequenced|pct low-diversity saturated|
|---|---|---:|---:|---:|---:|---:|---:|
|low_reads_50k|FALSE|188|87.9|483588|538|0.0|5.9|
|low_reads_50k|TRUE|26|12.1|25429|74|0.0|88.5|
|low_reads_100k|FALSE|182|85.0|500132|540|0.0|2.7|
|low_reads_100k|TRUE|32|15.0|36918|94|0.0|90.6|
|low_taxa_150|FALSE|186|86.9|485232|539|0.0|4.8|
|low_taxa_150|TRUE|28|13.1|27628|84|0.0|89.3|
|low_taxa_250|FALSE|177|82.7|511586|541|0.0|0.0|
|low_taxa_250|TRUE|37|17.3|39371|108|0.0|91.9|
|low_reads_or_taxa|FALSE|177|82.7|511586|541|0.0|0.0|
|low_reads_or_taxa|TRUE|37|17.3|39371|108|0.0|91.9|
|low_reads_and_taxa|FALSE|182|85.0|500132|540|0.0|2.7|
|low_reads_and_taxa|TRUE|32|15.0|36918|94|0.0|90.6|

## Samples needing attention

|sample|core|age kyr|reads|observed taxa|slope per 10k|class|
|---|---|---:|---:|---:|---:|---|
|LV3003047160|ST8|55.89|6403|24|0.00|insufficient_for_call|
|LV3003047074|ST8|100.67|8155|20|0.00|insufficient_for_call|
|LV3003047073|ST8|101.69|4586|12|0.00|insufficient_for_call|
|LV7008886741|GeoB25202_R1|59.00|36360|90|0.00|low_diversity_but_saturated|
|LV7008886705|GeoB25202_R1|59.74|82777|207|0.00|low_diversity_but_saturated|
|LV7008886781|GeoB25202_R1|120.02|24863|66|0.00|low_diversity_but_saturated|
|LV7008886757|GeoB25202_R1|123.12|41240|97|0.00|low_diversity_but_saturated|
|LV7008886228|GeoB25202_R2|120.02|41750|85|0.00|low_diversity_but_saturated|
|LV7008886769|GeoB25202_R2|123.12|112804|187|0.00|low_diversity_but_saturated|
|LV3003061473|ST13|117.09|119389|164|0.00|low_diversity_but_saturated|
|LV3003061568|ST13|121.64|22478|53|0.00|low_diversity_but_saturated|
|LV3003061555|ST13|126.44|37477|99|0.00|low_diversity_but_saturated|
|LV3003061527|ST13|133.78|14647|49|0.00|low_diversity_but_saturated|
|LV3003047028|ST8|4.75|11886|38|0.00|low_diversity_but_saturated|
|LV3003047141|ST8|39.15|6891|27|0.00|low_diversity_but_saturated|
|LV3003047140|ST8|46.09|29261|108|0.00|low_diversity_but_saturated|
|LV3003047127|ST8|47.62|96803|240|0.00|low_diversity_but_saturated|
|LV3003047161|ST8|54.98|7628|30|0.00|low_diversity_but_saturated|
|LV3003047125|ST8|72.87|107697|232|0.00|low_diversity_but_saturated|
|LV3003047124|ST8|74.23|69484|145|0.00|low_diversity_but_saturated|
|LV3003047123|ST8|75.59|25995|82|0.00|low_diversity_but_saturated|
|LV3003047122|ST8|76.95|43633|126|0.00|low_diversity_but_saturated|
|LV3003047105|ST8|82.45|62015|114|0.00|low_diversity_but_saturated|
|LV3003047104|ST8|83.84|88376|138|0.00|low_diversity_but_saturated|
|LV3003047090|ST8|90.76|110264|161|0.00|low_diversity_but_saturated|
|LV3003047089|ST8|92.14|141951|206|0.00|low_diversity_but_saturated|
|LV3003047088|ST8|93.52|46002|112|0.00|low_diversity_but_saturated|
|LV3003047087|ST8|94.90|45812|142|0.00|low_diversity_but_saturated|
|LV3003047064|ST8|96.28|21720|35|0.00|low_diversity_but_saturated|
|LV3003047078|ST8|98.76|13396|44|0.00|low_diversity_but_saturated|
|LV3003047077|ST8|99.71|46091|117|0.00|low_diversity_but_saturated|
|LV3003047162|ST8|102.70|24615|53|0.00|low_diversity_but_saturated|
|LV3003047063|ST8|103.69|24610|57|0.00|low_diversity_but_saturated|
|LV3003047062|ST8|104.69|39336|108|0.00|low_diversity_but_saturated|
|LV3003047061|ST8|105.69|39371|109|0.00|low_diversity_but_saturated|
|LV3003047060|ST8|106.69|46974|162|0.00|low_diversity_but_saturated|
|LV3003047059|ST8|107.68|79512|169|0.00|low_diversity_but_saturated|
## CLR/log transform note

CLR and log transforms do not recover unobserved taxa. They can reduce scale effects among observed taxa, but if a sample lacks reads or detected taxa, a transform cannot infer the missing biological diversity. That is why rarefaction and depth QC must be evaluated before deciding whether a transformed matrix is trustworthy for network construction.

## Outputs

- Sample metrics: `InputQC/rarefaction_depth_qc/results/tables/sample_rarefaction_depth_metrics.tsv`
- Rarefaction checkpoints: `InputQC/rarefaction_depth_qc/results/tables/rarefaction_checkpoint_table.tsv`
- Threshold summary: `InputQC/rarefaction_depth_qc/results/tables/threshold_depth_summary.tsv`
- QC class summary: `InputQC/rarefaction_depth_qc/results/tables/qc_depth_class_summary.tsv`
- Figures: `InputQC/rarefaction_depth_qc/results/figures/`
