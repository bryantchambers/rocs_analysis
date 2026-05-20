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
METADMG="/projects/caeg/apps/metaDMG_300326/metaDMG-cpp/metaDMG-cpp"

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

    # outbam="$EUK_OUT/${sample}.comp.filtered.bam"
    # outstat="$EUK_OUT/${sample}.comp.filtered.unicorn.refstats"
    # logfile="$LOGS/${sample}__unicorn_refstats.log"

    if [[ -s "$outbam" && -s "$outstat" ]]; then
        echo "Skipping $sample: outputs already exist"
        continue
    fi

    echo "Processing $sample"

    # "$UNICORN" refstats \
    #     -b "$bam" \
    #     -t "$THREADS" \
    #     -minreads 3 \
    #     --outbam "$outbam" \
    #     --outstat "$outstat" \
    #     --names "$TAX_PATH_NCBI/names.dmp" \
    #     --nodes "$TAX_PATH_NCBI/nodes.dmp" \
    #     > "$logfile" 2>&1
    
    log_step "Running unicorn alnfilter..."
    # parallel -j "$THREADSP" --colsep '\t' \
      "$UNICORN" alnfilt \
         -b "$bam" \
         -t $THREADS --mode $ALNFILT_MODE \
          --outbam $EUK_DIR/${sample}.alnfilt.bam \
          --outstat $EUK_DIR/${sample}.alnfilt.refstats \
          --minani $ALNFILT_MINANI --maxani $ALNFILT_MAXANI" \
      :::: "$SAMPLE_LIST"
  check_success "Unicorn alnfilt filtering"


  log_step "Running unicorn refstats filter..."
  parallel -j "$THREADSP" --colsep '\t' \
      "$UNICORN refstats \
         -b $EUK_DIR/${sample}.alnfilt.bam \
          -t $THREADS --minreads $UNICORN_MINREADS --outbam $EUK_DIR/${sample}.filtered.bam \
          --outstat $EUK_DIR/${sample}.filtered.refstats \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP" \
      :::: "$SAMPLE_LIST"
  check_success "Unicorn refstats final filtering"

  parallel -j "$THREADSP" --colsep '\t' '
    out_bam="'"$EUK_DIR"'/{1}.filtered.bam"
    out_stat="'"$EUK_DIR"'/{1}.filtered.refstats"
    [[ -s "$out_bam" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn BAM missing or empty: $out_bam" >> "'"$MISSING_DATA_LOG"'"
    [[ -s "$out_stat" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn refstats missing or empty: $out_stat" >> "'"$MISSING_DATA_LOG"'"
  ' :::: "$SAMPLE_LIST"

  log_step "Sorting merged BAM file for metaDMG..."
  parallel -j "$THREADSP" --colsep '\t' \
      "samtools sort -n -@ $THREADS -m 10G -o $EUK_DIR/{1}.sort.filtered.bam $EUK_DIR/{1}.filtered.bam" \
      :::: "$SAMPLE_LIST"
  check_success "Sorting BAM file"

  log_step "Running taxonomic classification with metaDMG..."
  parallel -j "$THREADSP" --colsep '\t' \
      "$METADMG lca \
        --names $TAX_NAMES_DMP \
        --nodes $TAX_NODES_DMP \
        --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
        --sim_score_low $METADMG_SIM_SCORE_LOW \
        --sim_score_high $METADMG_SIM_SCORE_HIGH \
        --how_many $METADMG_HOW_MANY \
        --weight_type $METADMG_WEIGHT_TYPE \
        --lca_rank $METADMG_LCA_RANK \
        --fix_ncbi $METADMG_FIX_NCBI \
        --threads $METADMG_THREADS \
        --filtered_acc2tax $EUK_DIR/{1}.acc2tax \
        --bam $EUK_DIR/{1}.sort.filtered.bam --out_prefix $EUK_DIR/{1}.sort.filtered" \
      :::: "$SAMPLE_LIST"
  check_success "Taxonomic classification"

  log_step "Running damage estimation with metaDMG..."
  parallel -j "$THREADSP" --colsep '\t' \
      "$METADMG dfit \
        $EUK_DIR/{1}.sort.filtered.bdamage.gz \
        --threads $METADMG_THREADS \
        --names $TAX_NAMES_DMP \
        --nodes $TAX_NODES_DMP \
        --showfits $METADMG_SHOWFITS \
        --nopt $METADMG_NOPT \
        --nbootstrap $METADMG_NBOOTSTRAP \
        --doboot $METADMG_DOBOOT \
        --seed $METADMG_SEED \
        --lib $METADMG_LIB \
        --out_prefix $EUK_DIR/{1}.sort.filtered" \
      :::: "$SAMPLE_LIST"
  check_success "Damage calculations done"

  log_step "Aggregating lca and dfit metaDMG..."
  parallel -j "$THREADSP" --colsep '\t' \
      "$METADMG aggregate \
        $EUK_DIR/{1}.sort.filtered.bdamage.gz \
        --names $TAX_NAMES_DMP \
        --nodes $TAX_NODES_DMP \
        --lcastat $EUK_DIR/{1}.sort.filtered.stat.gz \
        --dfit $EUK_DIR/{1}.sort.filtered.dfit.gz \
        --out_prefix $EUK_DIR/{1}.sort.filtered.agg" \
      :::: "$SAMPLE_LIST"
  check_success "Aggregation done"

  parallel -j "$THREADSP" --colsep '\t' '
    [[ -s "'"$EUK_DIR"'/{1}.sort.filtered.stat.gz" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=metaDMG stat.gz missing or empty" >> "'"$MISSING_DATA_LOG"'"
    [[ -s "'"$EUK_DIR"'/{1}.sort.filtered.dfit.gz" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=metaDMG dfit.gz missing or empty" >> "'"$MISSING_DATA_LOG"'"
    [[ -s "'"$EUK_DIR"'/{1}.sort.filtered.agg.stat.gz" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=metaDMG aggregate stat.gz missing or empty" >> "'"$MISSING_DATA_LOG"'"
  ' :::: "$SAMPLE_LIST"
  
  log_step "Unicorn per taxID statistics (family level)..."
    parallel -j "$THREADSP" --colsep '\t' \
      "$UNICORN taxstats \
         -b $EUK_DIR/{1}.filtered.bam \
          -t $THREADS \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP \
          --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
          -k $TAXSTATS_K --outstat $EUK_DIR/{1}.filtered.family.taxstats --minreads $TAXSTATS_MINREADS --rank $TAXSTATS_RANK_FAMILY" \
      :::: "$SAMPLE_LIST"
    check_success "Unicorn taxstats (family)"

    parallel -j "$THREADSP" --colsep '\t' '
      [[ -s "'"$EUK_DIR"'/{1}.filtered.family.taxstats" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn family taxstats missing or empty" >> "'"$MISSING_DATA_LOG"'"
    ' :::: "$SAMPLE_LIST"

    log_step "Unicorn per taxID statistics (genus level)..."
    parallel -j "$THREADSP" --colsep '\t' \
      "$UNICORN taxstats \
         -b $EUK_DIR/{1}.filtered.bam \
          -t $THREADS \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP \
          --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
          -k $TAXSTATS_K --outstat $EUK_DIR/{1}.filtered.genus.taxstats --minreads $TAXSTATS_MINREADS --rank $TAXSTATS_RANK_GENUS" \
      :::: "$SAMPLE_LIST"
    check_success "Unicorn taxstats (genus)"

    parallel -j "$THREADSP" --colsep '\t' '
      [[ -s "'"$EUK_DIR"'/{1}.filtered.genus.taxstats" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn genus taxstats missing or empty" >> "'"$MISSING_DATA_LOG"'"
    ' :::: "$SAMPLE_LIST"

    log_step "Unicorn per taxID statistics (species level)..."
    parallel -j "$THREADSP" --colsep '\t' \
      "$UNICORN taxstats \
         -b $EUK_DIR/{1}.filtered.bam \
          -t $THREADS \
          --names $TAX_NAMES_DMP \
          --nodes $TAX_NODES_DMP \
          --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz $ACC2TAX_EXTRA) \
          -k $TAXSTATS_K --outstat $EUK_DIR/{1}.filtered.spec.taxstats --minreads $TAXSTATS_MINREADS --rank $TAXSTATS_RANK_SPECIES" \
      :::: "$SAMPLE_LIST"
    check_success "Unicorn taxstats (species)"

    parallel -j "$THREADSP" --colsep '\t' '
      [[ -s "'"$EUK_DIR"'/{1}.filtered.spec.taxstats" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn species taxstats missing or empty" >> "'"$MISSING_DATA_LOG"'"
    ' :::: "$SAMPLE_LIST"



done
