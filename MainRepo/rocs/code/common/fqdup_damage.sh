#!/usr/bin/env bash
set -euo pipefail

# # Paths
# FQDUP="/maps/projects/caeg/people/ngm902/apps/repos/fqdup/build/fqdup"
# METADATA="/projects/caeg/people/ngm902/apps/repos/rocs/data/metadata_v5.tsv"
# READS_DIR="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/microbial/read-renamed"

# THREADS=48

# # Saltar header y iterar sobre labels
# tail -n +2 "$METADATA" | cut -f1 | while read -r label; do
    
#     fq="${READS_DIR}/${label}.fq.gz"
#     out_tsv="${label}_damage_profile.tsv"
#     out_json="${label}_damage.json"

#     # Check input exists
#     if [[ ! -f "$fq" ]]; then
#         echo "⚠️  Missing FASTQ for $label → skipping"
#         continue
#     fi

#     echo "▶ Running fqdup damage for $label"

#     "$FQDUP" damage \
#         -i "$fq" \
#         --library-type ds \
#         -p "$THREADS" \
#         --tsv "$out_tsv" \
#         --json "$out_json"

# done

# echo "DONE ✅"




INPUT_DIR="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/fqdup"
OUT="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/fqdup/fqdup_damage_summary.tsv"

echo "Collapsing fqdup damage profiles from $INPUT_DIR into $OUT"


printf '%s\n' \
"label	input	n_reads	library_type	library_type_auto	library_type_rescued	damage_status	d_max_5prime	d_max_3prime	d_max_combined	d_metamatch	damage_source	lambda_5prime	lambda_3prime	bg_5prime	bg_3prime	validated	artifact	cpg_ratio	dmax_ct5_cpg	dmax_ct5_noncpg	cpg_contrast	dipyr_contrast	heterogeneity_p	heterogeneity_detected	D	s_gt	s_oxog	channel_c_detected	short_asym_log2oe	short_log2oe	short_z	depurination_detected	depurination_enrichment_5prime	depurination_enrichment_3prime	depurination_rate_interior	bic_bias	bic_ds	bic_ss	bic_ct5_amp	bic_ga3_amp	bic_ga0_amp	bic_ct3_amp" \
> "$OUT"

shopt -s nullglob
for f in "$INPUT_DIR"/*.json; do
  jq -r '
    [
      (.input | split("/")[-1] | sub("\\.fq\\.gz$"; "")),
      .input,
      .n_reads,
      .library_type,
      .library_type_auto,
      .library_type_rescued,
      .damage_status,

      .deamination.d_max_5prime,
      .deamination.d_max_3prime,
      .deamination.d_max_combined,
      .deamination.d_metamatch,
      .deamination.source,
      .deamination.lambda_5prime,
      .deamination.lambda_3prime,
      .deamination.bg_5prime,
      .deamination.bg_3prime,
      .deamination.validated,
      .deamination.artifact,

      .deamination.cpg_like.cpg_ratio,
      .deamination.cpg_like.dmax_ct5_cpg,
      .deamination.cpg_like.dmax_ct5_noncpg,

      .deamination.context_deamination.cpg_contrast,
      .deamination.context_deamination.dipyr_contrast,
      .deamination.context_deamination.heterogeneity_p,
      .deamination.context_deamination.heterogeneity_detected,

      .complement_asymmetry.D,
      .complement_asymmetry.s_gt,
      .complement_asymmetry.s_oxog,
      .complement_asymmetry.channel_c_detected,

      .interior_ct_cluster.short_asym_log2oe,
      .interior_ct_cluster.short_log2oe,
      .interior_ct_cluster.short_z,

      .depurination.detected,
      .depurination.enrichment_5prime,
      .depurination.enrichment_3prime,
      .depurination.rate_interior,

      .bic.bias,
      .bic.ds,
      .bic.ss,
      .bic.ct5_amp,
      .bic.ga3_amp,
      .bic.ga0_amp,
      .bic.ct3_amp
    ] | @tsv
  ' "$f" >> "$OUT"
done

echo "DONE -> $OUT"