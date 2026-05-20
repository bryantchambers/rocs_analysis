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
# THREADS=10
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


    # if [[ -s "$EUK_DIR/${sample}.alnfilt.bam" && -s "$EUK_DIR/${sample}.alnfilt.refstats" ]]; then
    #     log_step "Skipping unicorn alnfilt for $sample: outputs already exist"
    # else
    #     if ! "$UNICORN" alnfilt \
    #       -b "$bam" \
    #       -t "$THREADS" --mode "$ALNFILT_MODE" \
    #       --outbam "$EUK_DIR/${sample}.alnfilt.bam" \
    #       --outstat "$EUK_DIR/${sample}.alnfilt.refstats" \
    #       --minani "$ALNFILT_MINANI" --maxani "$ALNFILT_MAXANI"; then

    #         log_step "WARNING: unicorn alnfilt failed for $sample. Skipping sample."
    #         echo -e "$sample\talnfilt_failed\t$bam" >> "$MISSING_DATA_LOG"
    #         continue
    #     fi
    # fi

    # # Optional safety check: if alnfilt BAM is missing/empty, skip
    # if [[ ! -s "$EUK_DIR/${sample}.alnfilt.bam" ]]; then
    #     log_step "WARNING: empty or missing alnfilt BAM for $sample. Skipping sample."
    #     echo -e "$sample\tempty_alnfilt_bam\t$bam" >> "$MISSING_DATA_LOG"
    #     continue
    # fi

    # if [[ -s "$EUK_DIR/${sample}.filtered.bam" && -s "$EUK_DIR/${sample}.filtered.refstats" ]]; then
    #     log_step "Skipping unicorn refstats for $sample: outputs already exist"
    # else
    #     if ! "$UNICORN" refstats \
    #       -b "$EUK_DIR/${sample}.alnfilt.bam" \
    #       -t "$THREADS" \
    #       --minreads "$UNICORN_MINREADS" \
    #       --outbam "$EUK_DIR/${sample}.filtered.bam" \
    #       --outstat "$EUK_DIR/${sample}.filtered.refstats" \
    #       --names "$TAX_NAMES_DMP" \
    #       --nodes "$TAX_NODES_DMP"; then

    #         log_step "WARNING: unicorn refstats failed for $sample. Skipping sample."
    #         echo -e "$sample\trefstats_failed\t$EUK_DIR/${sample}.alnfilt.bam" >> "$MISSING_DATA_LOG"
    #         continue
    #     fi
    # fi

    # # Safety check before taxstats
    # if [[ ! -s "$EUK_DIR/${sample}.filtered.bam" ]]; then
    #     log_step "WARNING: empty or missing filtered BAM for $sample. Skipping taxstats."
    #     echo -e "$sample\tempty_filtered_bam\t$EUK_DIR/${sample}.filtered.bam" >> "$MISSING_DATA_LOG"
    #     continue
    # fi




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



    # # ---------------------------
    # # TAXSTATS FAMILY
    # # ---------------------------
    # log_step "Unicorn taxstats - family for $sample"

    # if [[ -s "$EUK_DIR/${sample}.filtered.family.taxstats" ]]; then
    #     log_step "Skipping taxstats family for $sample: output exists"
    # else
    #     if ! "$UNICORN" taxstats \
    #       -b "$EUK_DIR/${sample}.filtered.bam" \
    #       -t "$THREADS" \
    #       --names "$TAX_NAMES_DMP" \
    #       --nodes "$TAX_NODES_DMP" \
    #       --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
    #       -k "$TAXSTATS_K" \
    #       --outstat "$EUK_DIR/${sample}.filtered.family.taxstats" \
    #       --minreads "$TAXSTATS_MINREADS" \
    #       --rank "$TAXSTATS_RANK_FAMILY"; then

    #         log_step "WARNING: taxstats family failed for $sample"
    #         echo -e "$sample\ttaxstats_family_failed" >> "$MISSING_DATA_LOG"
    #     fi
    # fi


    # # ---------------------------
    # # TAXSTATS GENUS
    # # ---------------------------
    # log_step "Unicorn taxstats - genus for $sample"

    # if [[ -s "$EUK_DIR/${sample}.filtered.genus.taxstats" ]]; then
    #     log_step "Skipping taxstats genus for $sample: output exists"
    # else
    #     if ! "$UNICORN" taxstats \
    #       -b "$EUK_DIR/${sample}.filtered.bam" \
    #       -t "$THREADS" \
    #       --names "$TAX_NAMES_DMP" \
    #       --nodes "$TAX_NODES_DMP" \
    #       --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
    #       -k "$TAXSTATS_K" \
    #       --outstat "$EUK_DIR/${sample}.filtered.genus.taxstats" \
    #       --minreads "$TAXSTATS_MINREADS" \
    #       --rank "$TAXSTATS_RANK_GENUS"; then

    #         log_step "WARNING: taxstats genus failed for $sample"
    #         echo -e "$sample\ttaxstats_genus_failed" >> "$MISSING_DATA_LOG"
    #     fi
    # fi


    # # ---------------------------
    # # TAXSTATS SPECIES
    # # ---------------------------
    # log_step "Unicorn taxstats - species for $sample"

    # if [[ -s "$EUK_DIR/${sample}.filtered.spec.taxstats" ]]; then
    #     log_step "Skipping taxstats species for $sample: output exists"
    # else
    #     if ! "$UNICORN" taxstats \
    #       -b "$EUK_DIR/${sample}.filtered.bam" \
    #       -t "$THREADS" \
    #       --names "$TAX_NAMES_DMP" \
    #       --nodes "$TAX_NODES_DMP" \
    #       --acc2tax <(zcat "$TAX_PATH_NCBI"/*.acc2taxid.gz "$ACC2TAX_EXTRA") \
    #       -k "$TAXSTATS_K" \
    #       --outstat "$EUK_DIR/${sample}.filtered.spec.taxstats" \
    #       --minreads "$TAXSTATS_MINREADS" \
    #       --rank "$TAXSTATS_RANK_SPECIES"; then

    #         log_step "WARNING: taxstats species failed for $sample"
    #         echo -e "$sample\ttaxstats_species_failed" >> "$MISSING_DATA_LOG"
    #     fi
    # fi

    
    
    
    log_step "Sorting merged BAM file for metaDMG..."

    filtered_bam="$EUK_DIR/${sample}.filtered.bam"
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

    if [[ ! -s "$filtered_bam" ]]; then
        log_step "WARNING: missing or empty filtered BAM for $sample. Skipping metaDMG."
        echo -e "$sample\tmissing_filtered_bam\t$filtered_bam" >> "$MISSING_DATA_LOG"
        continue
    fi

    if [[ -s "$sort_bam" ]]; then
        log_step "Skipping samtools sort for $sample: output already exists"
    else
        if ! samtools sort -n \
            -@ "$THREADS" \
            -m 10G \
            -o "$sort_bam" \
            "$filtered_bam"; then

            log_step "WARNING: samtools sort failed for $sample. Skipping metaDMG."
            echo -e "$sample\tsamtools_sort_failed\t$filtered_bam" >> "$MISSING_DATA_LOG"
            continue
        fi
    fi

    if [[ ! -s "$sort_bam" ]]; then
        log_step "WARNING: missing or empty sorted BAM for $sample. Skipping metaDMG."
        echo -e "$sample\tmissing_sort_bam\t$sort_bam" >> "$MISSING_DATA_LOG"
        continue
    fi

    log_step "Running taxonomic classification with metaDMG..."

    if [[ -s "$lca_stat" && -s "$lca_bdamage" && -s "$lca_lca" && -s "$lca_rlens" && -s "$lca_usedreads" && -s "$acc2tax" ]]; then
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
            continue
        fi
    fi

    if [[ ! -s "$lca_bdamage" || ! -s "$lca_stat" ]]; then
        log_step "WARNING: missing metaDMG lca outputs for $sample. Skipping dfit/aggregate."
        echo -e "$sample\tmissing_lca_outputs\t$EUK_DIR/${sample}.sort.filtered" >> "$MISSING_DATA_LOG"
        continue
    fi

    log_step "Running damage estimation with metaDMG..."

    if [[ -s "$dfit" && -s "$boot_stat" ]]; then
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
            continue
        fi
    fi

    if [[ ! -s "$dfit" ]]; then
        log_step "WARNING: missing dfit output for $sample. Skipping aggregate."
        echo -e "$sample\tmissing_dfit\t$dfit" >> "$MISSING_DATA_LOG"
        continue
    fi

    log_step "Damage calculations done"

    log_step "Aggregating lca and dfit metaDMG..."

    if [[ -s "$agg_stat" ]]; then
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
            continue
        fi
    fi

    log_step "Aggregation done"







done
