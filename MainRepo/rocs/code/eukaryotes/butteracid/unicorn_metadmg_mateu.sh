#!/usr/bin/env bash
set -euo pipefail

# Update the symbolic links to the BAM files if they don't already exist
shopt -s nullglob

# for f in /projects/caeg/people/ngm902/apps/repos/aeDNA/results/aligns/merge/*; do
#   dest="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/bam/$(basename "$f")"
  
#   if [[ -e "$dest" ]]; then
#     echo "Skipping $(basename "$f") (already exists)"
#     continue
#   fi
  
#   ln -s "$f" "$dest"
# done


for f in /projects/caeg/people/ngm902/apps/aeDNA/results/aligns/merge/*; do
  dest="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/bam/$(basename "$f")"
  
  if [[ -e "$dest" ]]; then
    echo "Skipping $(basename "$f") (already exists)"
    continue
  fi
  
  ln -s "$f" "$dest"
done





# ---------------------------
# Paths and executables
# ---------------------------
EUK_IN="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/bam"
EUK_OUT="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/unicorn_rerun"
EUK_DIR="$EUK_OUT"
LOGS="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/unicorn_rerun/logs"
SLURM_LOGS="$LOGS/slurm"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

UNICORN="/projects/caeg/apps/unicorn_300326/unicorn/unicorn"
METADMG="/projects/caeg/people/ngm902/apps/repos/metaDMG-cpp/metaDMG-cpp"

DB_PATH="/datasets/caeg_dataset/references/ncbi/20250530/data/wgs_eukaryota"
DB_PATH_CLEAN="/datasets/caeg_dataset/references/ncbi/20250530/data/"
TAX_PATH_NCBI="/datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi"
TAXDUMP_DIR="$TAX_PATH_NCBI/taxdump"
TAX_NAMES_DMP="$TAXDUMP_DIR/names.dmp"
TAX_NODES_DMP="$TAXDUMP_DIR/nodes.dmp"
ACC2TAX_EXTRA="/datasets/caeg_dataset/references/phylo_norway/phylo_norway.acc2taxid.gz"

# ---------------------------
# Runtime controls
# ---------------------------
THREADS="${THREADS:-10}"
THREADSP="$THREADS"
SAMPLE_LIST="/projects/caeg/people/ngm902/apps/repos/rocs/code/eukaryotes/butteracid/sample_list.tsv"
MISSING_DATA_LOG="$LOGS/missing_data.log"

# ---------------------------
# SLURM controls
# ---------------------------
MAX_PARALLEL_JOBS="${MAX_PARALLEL_JOBS:-60}"
SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-$THREADS}"
SLURM_MEM="${SLURM_MEM:-150G}"
SLURM_PARTITION="${SLURM_PARTITION:-comppriority}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-prio}"
SLURM_TIME="${SLURM_TIME:-24:00:00}"
SBATCH_EXTRA_ARGS="${SBATCH_EXTRA_ARGS:-}"

# ---------------------------
# Unicorn / metaDMG parameters
# ---------------------------
ALNFILT_MODE="ALLTOP"
ALNFILT_MINANI="93.00"
ALNFILT_MAXANI="100.00"
UNICORN_MINREADS="3"

METADMG_SIM_SCORE_LOW="0.93"
METADMG_SIM_SCORE_HIGH="1.0"
METADMG_HOW_MANY="15"
METADMG_WEIGHT_TYPE="1"
METADMG_LCA_RANK="family"
METADMG_FIX_NCBI="0"
METADMG_THREADS="$THREADS"

METADMG_SHOWFITS="2"
METADMG_NOPT="10"
METADMG_NBOOTSTRAP="20"
METADMG_DOBOOT="1"
METADMG_SEED="1234"
METADMG_LIB="ds"

TAXSTATS_K="17"
TAXSTATS_MINREADS="$UNICORN_MINREADS"
TAXSTATS_RANK_FAMILY="family"
TAXSTATS_RANK_GENUS="genus"
TAXSTATS_RANK_SPECIES="species"

mkdir -p "$EUK_OUT" "$LOGS" "$SLURM_LOGS"
shopt -s nullglob
ACTIVE_JOB_IDS=()

log_step() {
    echo "[$(date)] $1"
}

all_outputs_exist() {
    local path
    for path in "$@"; do
        if [[ ! -s "$path" ]]; then
            return 1
        fi
    done
    return 0
}

log_missing_outputs() {
    local sample="$1"
    local step_name="$2"
    shift 2

    local path
    for path in "$@"; do
        if [[ ! -s "$path" ]]; then
            log_step "WARNING: missing output for $sample [$step_name]: $path"
            echo -e "$sample\t${step_name}_missing\t$path" >> "$MISSING_DATA_LOG"
        fi
    done
}

set_sample_paths() {
    local sample="$1"

    alnfilt_bam="$EUK_DIR/${sample}.alnfilt.bam"
    alnfilt_refstats="$EUK_DIR/${sample}.alnfilt.refstats"
    filtered_bam="$EUK_DIR/${sample}.filtered.bam"
    filtered_refstats="$EUK_DIR/${sample}.filtered.refstats"
    family_taxstats="$EUK_DIR/${sample}.filtered.family.taxstats"
    genus_taxstats="$EUK_DIR/${sample}.filtered.genus.taxstats"
    species_taxstats="$EUK_DIR/${sample}.filtered.spec.taxstats"
    sort_bam="$EUK_DIR/${sample}.sort.filtered.bam"
    acc2tax="$EUK_DIR/${sample}.acc2tax"
    lca_stat="$EUK_DIR/${sample}.sort.filtered.stat.gz"
    lca_bdamage="$EUK_DIR/${sample}.sort.filtered.bdamage.gz"
    lca_lca="$EUK_DIR/${sample}.sort.filtered.lca.gz"
    lca_rlens="$EUK_DIR/${sample}.sort.filtered.rlens.gz"
    lca_usedreads="$EUK_DIR/${sample}.sort.filtered.usedreads.bam"
    dfit="$EUK_DIR/${sample}.sort.filtered.dfit.gz"
    boot_stat="$EUK_DIR/${sample}.sort.filtered.boot.stat.gz"
    agg_stat="$EUK_DIR/${sample}.sort.filtered.agg.stat.gz"
}

sample_is_complete() {
    local sample="$1"
    set_sample_paths "$sample"

    all_outputs_exist \
        "$filtered_bam" \
        "$filtered_refstats" \
        "$family_taxstats" \
        "$genus_taxstats" \
        "$species_taxstats" \
        "$sort_bam" \
        "$acc2tax" \
        "$lca_stat" \
        "$lca_bdamage" \
        "$lca_lca" \
        "$lca_rlens" \
        "$lca_usedreads" \
        "$dfit" \
        "$boot_stat" \
        "$agg_stat"
}

process_sample() {
    local bam="$1"
    local sample
    sample=$(basename "$bam" .bam)
    set_sample_paths "$sample"

    log_step "Processing $sample"

    if all_outputs_exist "$alnfilt_bam"; then
        log_step "Skipping unicorn alnfilt for $sample: output already exists"
    else
        if ! "$UNICORN" alnfilt \
          -b "$bam" \
          -t "$THREADS" --mode "$ALNFILT_MODE" \
          --outbam "$alnfilt_bam" \
          --outstat "$alnfilt_refstats" \
          --minani "$ALNFILT_MINANI" --maxani "$ALNFILT_MAXANI"; then
            log_step "WARNING: unicorn alnfilt failed for $sample. Skipping sample."
            echo -e "$sample\talnfilt_failed\t$bam" >> "$MISSING_DATA_LOG"
            return 0
        fi
    fi

    if ! all_outputs_exist "$alnfilt_bam"; then
        log_missing_outputs "$sample" "alnfilt" "$alnfilt_bam"
        log_step "WARNING: incomplete unicorn alnfilt outputs for $sample. Skipping sample."
        return 0
    fi

    if all_outputs_exist "$filtered_bam" "$filtered_refstats"; then
        log_step "Skipping unicorn refstats for $sample: outputs already exist"
    else
        if ! "$UNICORN" refstats \
          -b "$alnfilt_bam" \
          -t "$THREADS" \
          --minreads "$UNICORN_MINREADS" \
          --outbam "$filtered_bam" \
          --outstat "$filtered_refstats" \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP"; then
            log_step "WARNING: unicorn refstats failed for $sample. Skipping sample."
            echo -e "$sample\trefstats_failed\t$alnfilt_bam" >> "$MISSING_DATA_LOG"
            return 0
        fi
    fi

    if ! all_outputs_exist "$filtered_bam" "$filtered_refstats"; then
        log_missing_outputs "$sample" "refstats" "$filtered_bam" "$filtered_refstats"
        log_step "WARNING: incomplete unicorn refstats outputs for $sample. Skipping sample."
        return 0
    fi

    log_step "Unicorn taxstats - family for $sample"
    if all_outputs_exist "$family_taxstats"; then
        log_step "Skipping taxstats family for $sample: output exists"
    else
        if ! "$UNICORN" taxstats \
          -b "$filtered_bam" \
          -t "$THREADS" \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP" \
          --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
          -k "$TAXSTATS_K" \
          --outstat "$family_taxstats" \
          --minreads "$TAXSTATS_MINREADS" \
          --rank "$TAXSTATS_RANK_FAMILY"; then
            log_step "WARNING: taxstats family failed for $sample"
            echo -e "$sample\ttaxstats_family_failed\t$family_taxstats" >> "$MISSING_DATA_LOG"
        fi
    fi
    if ! all_outputs_exist "$family_taxstats"; then
        log_missing_outputs "$sample" "taxstats_family" "$family_taxstats"
    fi

    log_step "Unicorn taxstats - genus for $sample"
    if all_outputs_exist "$genus_taxstats"; then
        log_step "Skipping taxstats genus for $sample: output exists"
    else
        if ! "$UNICORN" taxstats \
          -b "$filtered_bam" \
          -t "$THREADS" \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP" \
          --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
          -k "$TAXSTATS_K" \
          --outstat "$genus_taxstats" \
          --minreads "$TAXSTATS_MINREADS" \
          --rank "$TAXSTATS_RANK_GENUS"; then
            log_step "WARNING: taxstats genus failed for $sample"
            echo -e "$sample\ttaxstats_genus_failed\t$genus_taxstats" >> "$MISSING_DATA_LOG"
        fi
    fi
    if ! all_outputs_exist "$genus_taxstats"; then
        log_missing_outputs "$sample" "taxstats_genus" "$genus_taxstats"
    fi

    log_step "Unicorn taxstats - species for $sample"
    if all_outputs_exist "$species_taxstats"; then
        log_step "Skipping taxstats species for $sample: output exists"
    else
        if ! "$UNICORN" taxstats \
          -b "$filtered_bam" \
          -t "$THREADS" \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP" \
          --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
          -k "$TAXSTATS_K" \
          --outstat "$species_taxstats" \
          --minreads "$TAXSTATS_MINREADS" \
          --rank "$TAXSTATS_RANK_SPECIES"; then
            log_step "WARNING: taxstats species failed for $sample"
            echo -e "$sample\ttaxstats_species_failed\t$species_taxstats" >> "$MISSING_DATA_LOG"
        fi
    fi
    if ! all_outputs_exist "$species_taxstats"; then
        log_missing_outputs "$sample" "taxstats_species" "$species_taxstats"
    fi

    log_step "Sorting merged BAM file for metaDMG..."
    if all_outputs_exist "$sort_bam"; then
        log_step "Skipping samtools sort for $sample: output already exists"
    else
        if ! samtools sort -n \
          -@ "$THREADS" \
          -m 10G \
          -o "$sort_bam" \
          "$filtered_bam"; then
            log_step "WARNING: samtools sort failed for $sample. Skipping metaDMG."
            echo -e "$sample\tsamtools_sort_failed\t$filtered_bam" >> "$MISSING_DATA_LOG"
            return 0
        fi
    fi

    if ! all_outputs_exist "$sort_bam"; then
        log_missing_outputs "$sample" "samtools_sort" "$sort_bam"
        log_step "WARNING: incomplete samtools sort outputs for $sample. Skipping metaDMG."
        return 0
    fi

    log_step "Running taxonomic classification with metaDMG..."
    if all_outputs_exist "$acc2tax" "$lca_stat" "$lca_bdamage" "$lca_lca" "$lca_rlens" "$lca_usedreads"; then
        log_step "Skipping metaDMG lca for $sample: outputs already exist"
    else
        if ! "$METADMG" lca \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP" \
          --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
          --sim_score_low "$METADMG_SIM_SCORE_LOW" \
          --sim_score_high "$METADMG_SIM_SCORE_HIGH" \
          --how_many "$METADMG_HOW_MANY" \
          --weight_type "$METADMG_WEIGHT_TYPE" \
          --lca_rank "$METADMG_LCA_RANK" \
          --fix_ncbi "$METADMG_FIX_NCBI" \
          --threads "$METADMG_THREADS" \
          --filtered_acc2tax "$acc2tax" \
          --bam "$sort_bam" \
          --out_prefix "$EUK_DIR/${sample}.sort.filtered"; then
            log_step "WARNING: metaDMG lca failed for $sample. Skipping dfit/aggregate."
            echo -e "$sample\tmetadmg_lca_failed\t$sort_bam" >> "$MISSING_DATA_LOG"
            return 0
        fi
    fi

    if ! all_outputs_exist "$acc2tax" "$lca_stat" "$lca_bdamage" "$lca_lca" "$lca_rlens" "$lca_usedreads"; then
        log_missing_outputs "$sample" "metadmg_lca" "$acc2tax" "$lca_stat" "$lca_bdamage" "$lca_lca" "$lca_rlens" "$lca_usedreads"
        log_step "WARNING: incomplete metaDMG lca outputs for $sample. Skipping dfit/aggregate."
        return 0
    fi

    log_step "Running damage estimation with metaDMG..."
    if all_outputs_exist "$dfit" "$boot_stat"; then
        log_step "Skipping metaDMG dfit for $sample: outputs already exist"
    else
        if ! "$METADMG" dfit \
          "$lca_bdamage" \
          --threads "$METADMG_THREADS" \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP" \
          --showfits "$METADMG_SHOWFITS" \
          --nopt "$METADMG_NOPT" \
          --nbootstrap "$METADMG_NBOOTSTRAP" \
          --doboot "$METADMG_DOBOOT" \
          --seed "$METADMG_SEED" \
          --lib "$METADMG_LIB" \
          --out_prefix "$EUK_DIR/${sample}.sort.filtered"; then
            log_step "WARNING: metaDMG dfit failed for $sample. Skipping aggregate."
            echo -e "$sample\tmetadmg_dfit_failed\t$lca_bdamage" >> "$MISSING_DATA_LOG"
            return 0
        fi
    fi

    if ! all_outputs_exist "$dfit" "$boot_stat"; then
        log_missing_outputs "$sample" "metadmg_dfit" "$dfit" "$boot_stat"
        log_step "WARNING: incomplete metaDMG dfit outputs for $sample. Skipping aggregate."
        return 0
    fi

    log_step "Damage calculations done"
    log_step "Aggregating lca and dfit metaDMG..."

    if all_outputs_exist "$agg_stat"; then
        log_step "Skipping metaDMG aggregate for $sample: output already exists"
    else
        if ! "$METADMG" aggregate \
          "$lca_bdamage" \
          --names "$TAX_NAMES_DMP" \
          --nodes "$TAX_NODES_DMP" \
          --lcastat "$lca_stat" \
          --dfit "$dfit" \
          --out_prefix "$EUK_DIR/${sample}.sort.filtered.agg"; then
            log_step "WARNING: metaDMG aggregate failed for $sample."
            echo -e "$sample\tmetadmg_aggregate_failed\t$lca_bdamage" >> "$MISSING_DATA_LOG"
            return 0
        fi
    fi

    if ! all_outputs_exist "$agg_stat"; then
        log_missing_outputs "$sample" "metadmg_aggregate" "$agg_stat"
        log_step "WARNING: incomplete metaDMG aggregate outputs for $sample."
        return 0
    fi

    log_step "Aggregation done for $sample"
}

refresh_active_job_ids() {
    local joined_ids

    if (( ${#ACTIVE_JOB_IDS[@]} == 0 )); then
        return 0
    fi

    joined_ids=$(IFS=,; echo "${ACTIVE_JOB_IDS[*]}")
    mapfile -t ACTIVE_JOB_IDS < <(squeue -a -h -j "$joined_ids" -o "%A")
}

wait_for_available_slot() {
    while true; do
        refresh_active_job_ids
        if (( ${#ACTIVE_JOB_IDS[@]} < MAX_PARALLEL_JOBS )); then
            break
        fi
        sleep 5
    done
}

build_sbatch_command() {
    local sample="$1"
    local bam="$2"
    local stdout_log="$SLURM_LOGS/job_unicorn_${sample}_%u_%j.out"
    local stderr_log="$SLURM_LOGS/job_unicorn_${sample}_%u_%j.err"

    SBATCH_CMD=(
        sbatch
        --parsable
        --job-name="unicorn_${sample}"
        --partition="$SLURM_PARTITION"
        --account="$SLURM_ACCOUNT"
        --cpus-per-task="$SLURM_CPUS_PER_TASK"
        --time="$SLURM_TIME"
        --mem="$SLURM_MEM"
        --output="$stdout_log"
        --error="$stderr_log"
        --wrap="bash \"$SCRIPT_DIR/$(basename "$0")\" --worker \"$bam\""
    )

    if [[ -n "$SBATCH_EXTRA_ARGS" ]]; then
        # shellcheck disable=SC2206
        extra_args=( $SBATCH_EXTRA_ARGS )
        SBATCH_CMD+=("${extra_args[@]}")
    fi
}

submit_all_samples() {
    local bam
    local submitted=0
    local skipped_complete=0
    local job_id

    for bam in "$EUK_IN"/*.bam; do
        local sample
        sample=$(basename "$bam" .bam)

        if sample_is_complete "$sample"; then
            log_step "Skipping submission for $sample: all expected outputs already exist"
            skipped_complete=$((skipped_complete + 1))
            continue
        fi

        wait_for_available_slot
        build_sbatch_command "$sample" "$bam"
        log_step "Submitting $sample via sbatch (max parallel jobs in queue: $MAX_PARALLEL_JOBS)"
        job_id=$("${SBATCH_CMD[@]}")
        job_id=${job_id%%;*}
        ACTIVE_JOB_IDS+=("$job_id")
        log_step "Submitted $sample as job $job_id"
        submitted=$((submitted + 1))
    done

    if (( submitted == 0 && skipped_complete == 0 )); then
        log_step "No BAM files found in $EUK_IN"
        return 0
    fi

    log_step "Submitted $submitted sample jobs via sbatch"
    log_step "Skipped $skipped_complete already-complete samples"
}

main() {
    if [[ "${1:-}" == "--worker" ]]; then
        if [[ $# -ne 2 ]]; then
            echo "Usage: $0 --worker /path/to/sample.bam" >&2
            exit 1
        fi
        process_sample "$2"
        return 0
    fi

    submit_all_samples
}

main "$@"
