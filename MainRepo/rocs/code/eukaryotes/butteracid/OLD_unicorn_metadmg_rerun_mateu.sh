#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# Paths and executables
# ---------------------------
EUK_IN="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/bam"
EUK_OUT="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/unicorn_rerun"
EUK_DIR="$EUK_OUT"
LOGS="/projects/caeg/people/ngm902/apps/repos/rocs/taxonomic-profiling/eukaryotes/unicorn_rerun/logs"

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
# Runtime controls and placeholders
# ---------------------------
THREADS=10
THREADSP="$THREADS"
# SAMPLE_LIST="$EUK_DIR/sample_list.tsv"
SAMPLE_LIST="/projects/caeg/people/ngm902/apps/repos/rocs/code/eukaryotes/butteracid/sample_list.tsv"
MISSING_DATA_LOG="$LOGS/missing_data.log"



# ---------------------------
# Logging helpers
# ---------------------------
log_step() {
    echo "[$(date)] $1"
}

check_success() {
    local step_name="$1"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] ${step_name} failed" >&2
        exit 1
    fi
}


# ---------------------------
# Step 5a: Unicorn alnfilt parameters
# ---------------------------
ALNFILT_MODE="ALLTOP"
ALNFILT_MINANI="93.00"
ALNFILT_MAXANI="100.00"

# ---------------------------
# Step 5a/b: Unicorn refstats parameters
# ---------------------------
UNICORN_MINREADS="3"

# ---------------------------
# Step 5a: metaDMG lca parameters
# ---------------------------
METADMG_SIM_SCORE_LOW="0.93"
METADMG_SIM_SCORE_HIGH="1.0"
METADMG_HOW_MANY="15"
METADMG_WEIGHT_TYPE="1"
METADMG_LCA_RANK="family"
METADMG_FIX_NCBI="0"
METADMG_THREADS="10"

# ---------------------------
# Step 5a: metaDMG dfit parameters
# ---------------------------
METADMG_SHOWFITS="2"
METADMG_NOPT="10"
METADMG_NBOOTSTRAP="20"
METADMG_DOBOOT="1"
METADMG_SEED="1234"
METADMG_LIB="ds"

# ---------------------------
# Step 5a: metaDMG aggregate parameters
# ---------------------------
# (no extra params; uses lca + dfit outputs with shared naming)

# ---------------------------
# Step 5a: Unicorn taxstats parameters
# ---------------------------
TAXSTATS_K="17"
TAXSTATS_MINREADS="$UNICORN_MINREADS"
TAXSTATS_RANK_FAMILY="family"
TAXSTATS_RANK_GENUS="genus"
TAXSTATS_RANK_SPECIES="species"



mkdir -p "$EUK_OUT" "$LOGS"

shopt -s nullglob

for bam in "$EUK_IN"/*.bam; do
    
    sample=$(basename "$bam" .bam)

    echo "Processing $sample"

    if [[ -s "$EUK_DIR/${sample}.alnfilt.bam" && -s "$EUK_DIR/${sample}.alnfilt.refstats" ]]; then
        log_step "Skipping unicorn alnfilt for $sample: outputs already exist"
    else
        "$UNICORN" alnfilt \
          -b "$bam" \
          -t $THREADS --mode $ALNFILT_MODE \
          --outbam $EUK_DIR/${sample}.alnfilt.bam \
          --outstat $EUK_DIR/${sample}.alnfilt.refstats \
          --minani $ALNFILT_MINANI --maxani $ALNFILT_MAXANI
    fi

    if [[ -s "$EUK_DIR/${sample}.filtered.bam" && -s "$EUK_DIR/${sample}.filtered.refstats" ]]; then
        log_step "Skipping unicorn refstats for $sample: outputs already exist"
    else
        "$UNICORN" refstats \
          -b $EUK_DIR/${sample}.alnfilt.bam \
          -t $THREADS --minreads $UNICORN_MINREADS --outbam $EUK_DIR/${sample}.filtered.bam \
          --outstat $EUK_DIR/${sample}.filtered.refstats \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP
    fi

    # log_step "Sorting merged BAM file for metaDMG..."
      
    # if [[ -s "$EUK_DIR/${sample}.sort.filtered.bam" ]]; then
    #     log_step "Skipping samtools sort for $sample: output already exists"
    # else
    #     samtools sort -n -@ $THREADS -m 10G -o $EUK_DIR/${sample}.sort.filtered.bam $EUK_DIR/${sample}.filtered.bam
    # fi
    
    # log_step "Running taxonomic classification with metaDMG..."

    # if [[ -s "$EUK_DIR/${sample}.sort.filtered.stat.gz" && -s "$EUK_DIR/${sample}.sort.filtered.bdamage.gz" && -s "$EUK_DIR/${sample}.acc2tax" ]]; then
    #     log_step "Skipping metaDMG lca for $sample: outputs already exist"
    # else
    #     "$METADMG" lca \
    #       --names $TAX_NAMES_DMP \
    #       --nodes $TAX_NODES_DMP \
    #       --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
    #       --sim_score_low $METADMG_SIM_SCORE_LOW \
    #       --sim_score_high $METADMG_SIM_SCORE_HIGH \
    #       --how_many $METADMG_HOW_MANY \
    #       --weight_type $METADMG_WEIGHT_TYPE \
    #       --lca_rank $METADMG_LCA_RANK \
    #       --fix_ncbi $METADMG_FIX_NCBI \
    #       --threads $METADMG_THREADS \
    #       --filtered_acc2tax $EUK_DIR/${sample}.acc2tax \
    #       --bam $EUK_DIR/${sample}.sort.filtered.bam --out_prefix $EUK_DIR/${sample}.sort.filtered
    # fi


    # log_step "Running damage estimation with metaDMG..."

    # if [[ -s "$EUK_DIR/${sample}.sort.filtered.dfit.gz" ]]; then
    #     log_step "Skipping metaDMG dfit for $sample: output already exists"
    # else
    #     "$METADMG" dfit \
    #       $EUK_DIR/${sample}.sort.filtered.bdamage.gz \
    #       --threads $METADMG_THREADS \
    #       --names $TAX_NAMES_DMP \
    #       --nodes $TAX_NODES_DMP \
    #       --showfits $METADMG_SHOWFITS \
    #       --nopt $METADMG_NOPT \
    #       --nbootstrap $METADMG_NBOOTSTRAP \
    #       --doboot $METADMG_DOBOOT \
    #       --seed $METADMG_SEED \
    #       --lib $METADMG_LIB \
    #       --out_prefix $EUK_DIR/${sample}.sort.filtered
    # fi
          
    # log_step "Damage calculations done"

    # log_step "Aggregating lca and dfit metaDMG..."

    # if [[ -s "$EUK_DIR/${sample}.sort.filtered.agg.stat.gz" ]]; then
    #     log_step "Skipping metaDMG aggregate for $sample: output already exists"
    # else
    #     "$METADMG" aggregate \
    #       $EUK_DIR/${sample}.sort.filtered.bdamage.gz \
    #       --names $TAX_NAMES_DMP \
    #       --nodes $TAX_NODES_DMP \
    #       --lcastat $EUK_DIR/${sample}.sort.filtered.stat.gz \
    #       --dfit $EUK_DIR/${sample}.sort.filtered.dfit.gz \
    #       --out_prefix $EUK_DIR/${sample}.sort.filtered.agg
    # fi

    # log_step "Aggregation done"





    log_step "Unicorn per taxID statistics - family level..."
        
    if [[ -s "$EUK_DIR/${sample}.filtered.family.taxstats" ]]; then
        log_step "Skipping unicorn taxstats family for $sample: output already exists"
    else
        "$UNICORN" taxstats \
          -b $EUK_DIR/${sample}.filtered.bam \
          -t $THREADS \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP \
          --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
          -k $TAXSTATS_K --outstat $EUK_DIR/${sample}.filtered.family.taxstats --minreads $TAXSTATS_MINREADS --rank $TAXSTATS_RANK_FAMILY
    fi

    log_step "Unicorn taxstats - family"

        
    log_step "Unicorn per taxID statistics - genus level..."

    if [[ -s "$EUK_DIR/${sample}.filtered.genus.taxstats" ]]; then
        log_step "Skipping unicorn taxstats genus for $sample: output already exists"
    else
        "$UNICORN" taxstats \
          -b $EUK_DIR/${sample}.filtered.bam \
          -t $THREADS \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP \
          --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
          -k $TAXSTATS_K --outstat $EUK_DIR/${sample}.filtered.genus.taxstats --minreads $TAXSTATS_MINREADS --rank $TAXSTATS_RANK_GENUS
    fi
          

    log_step "Unicorn per taxID statistics - species level..."

    if [[ -s "$EUK_DIR/${sample}.filtered.spec.taxstats" ]]; then
        log_step "Skipping unicorn taxstats species for $sample: output already exists"
    else
        "$UNICORN" taxstats \
          -b $EUK_DIR/${sample}.filtered.bam \
          -t $THREADS \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP \
          --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
          -k $TAXSTATS_K --outstat $EUK_DIR/${sample}.filtered.spec.taxstats --minreads $TAXSTATS_MINREADS --rank $TAXSTATS_RANK_SPECIES
    fi

    log_step "Unicorn taxstats - species"


done
