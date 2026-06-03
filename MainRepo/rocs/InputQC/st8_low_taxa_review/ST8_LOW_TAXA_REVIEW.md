# ST8 Low-Taxa Sample Review

- Generated: 2026-05-19 11:52:55 CEST
- Definition: ST8 samples with `detected_taxa < 250`
- Progress log: `InputQC/st8_low_taxa_review/results/st8_low_taxa_review.log`

## Main findings

There are `27` ST8 low-taxa samples under the `<250 detected taxa` threshold.

Their library concentration is substantially lower than the rest of ST8 and the overall dataset, but not uniformly low; several low-taxa samples still have moderate library concentration. That argues for a mixed process: DNA/library yield contributes strongly, but it probably does not explain every low-taxa ST8 sample by itself.

## DNA/library concentration comparison

|group|n|median library concentration|IQR|median reads|median detected taxa|median Shannon|median age kyr|
|---|---:|---:|---:|---:|---:|---:|---:|
|ST8 low taxa (<250)|27|36.23|24.88-60.30|39371|109|4.00|92.14|
|ST8 not low taxa|91|99.94|70.49-131.87|361920|508|5.39|33.28|
|non-ST8|100|93.77|46.02-185.32|619230|558|5.05|59.37|
|all samples|218|84.65|48.48-136.78|405218|497|5.20|50.52|

## Most recurrent taxa

|taxon|samples present|pct samples|total counts|domain|phylum|functional group|signal source|
|---|---:|---:|---:|---|---|---|---|
|S__GCA_011055585.1|27|100.0|85257|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__GCA_016838845.1|27|100.0|54034|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__3300024058_12|26|96.3|38374|d__Bacteria|p__Proteobacteria|Core_heterotrophy|pelagic_export|
|S__GCA_013911135.1|25|92.6|62542|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__GCA_013390475.1|25|92.6|19114|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__GCA_008974855.1|24|88.9|25683|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__3300006311_3|23|85.2|12766|d__Bacteria|p__Chloroflexota|Dehalococcoidia|diagenetic_in_situ|
|S__GCA_016838795.1|22|81.5|32592|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__3300026253_28|22|81.5|7741|d__Bacteria|p__Acidobacteriota|Core_heterotrophy|pelagic_export|
|S__3300024516_13|21|77.8|15696|d__Bacteria|p__Proteobacteria|Core_heterotrophy|pelagic_export|
|S__GCA_905182815.1|21|77.8|6936|d__Bacteria|p__Proteobacteria|Pelagic_heterotroph|pelagic_export|
|S__GCA_905182865.1|21|77.8|6560|d__Bacteria|p__SAR324|Core_heterotrophy|pelagic_export|
|S__GCA_011049715.1|20|74.1|19066|d__Bacteria|p__Desulfobacterota|Sulfate_reduction_diagenetic|diagenetic_in_situ|
|S__GCA_017368835.1|20|74.1|7452|d__Bacteria|p__Planctomycetota|Anammox|particle_associated|
|S__GCA_001437625.1|20|74.1|5212|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__GCA_014384505.1|20|74.1|5128|d__Bacteria|p__Bacteroidota|Particle_heterotroph|particle_associated|
|S__GCA_905181595.1|20|74.1|4771|d__Bacteria|p__Proteobacteria|Pelagic_heterotroph|pelagic_export|
|S__3300024516_38|19|70.4|7721|d__Bacteria|p__KSB1|Unknown|uncertain|
|S__3300024058_22|19|70.4|6234|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|
|S__GCF_002836475.1|19|70.4|5970|d__Bacteria|p__Bacteroidota|Core_heterotrophy|pelagic_export|
|S__3300027752_17|19|70.4|5900|d__Bacteria|p__Proteobacteria|Core_heterotrophy|pelagic_export|
|S__GCA_015659035.1|19|70.4|5632|d__Bacteria|p__Acidobacteriota|Core_heterotrophy|pelagic_export|
|S__3300026253_8|19|70.4|5294|d__Bacteria|p__Planctomycetota|Core_heterotrophy|pelagic_export|
|S__TARA_ARC_108_MAG_00146|19|70.4|5247|d__Bacteria|p__Verrucomicrobiota|Particle_heterotroph|particle_associated|
|S__GCA_013215285.1|18|66.7|10785|d__Archaea|p__Thermoproteota|Nitrification_AOA|pelagic_export|

## Read distribution across detected taxa

Across ST8 low-taxa samples, the median top taxon contains `11.7%` of reads, the median top 5 taxa contain `36.0%`, and the median top 10 taxa contain `49.7%`.

The median sample needs `11` taxa to explain 50% of reads, `45` taxa to explain 80%, and `75` taxa to explain 90%. Median Pielou evenness is `0.87`, so most samples are not collapsed into a single dominant taxon, although several sparse samples are strongly concentrated in their top few taxa.

## Outputs

- Sample table: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_sample_table.tsv`
- Taxa summary: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_taxa_summary.tsv`
- Top taxa per sample: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_top_taxa_per_sample.tsv`
- Function summary: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_function_summary.tsv`
- Read distribution: `InputQC/st8_low_taxa_review/results/tables/st8_low_taxa_read_distribution.tsv`
- Figures: `InputQC/st8_low_taxa_review/results/figures/`
