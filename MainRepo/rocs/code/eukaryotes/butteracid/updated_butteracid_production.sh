#!/bin/bash
set -euo pipefail

# ---------------------------
# Load configuration
# ---------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 CONFIG_FILE"
    exit 1
fi

CONFIG_FILE="$1"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "[ERROR] Config file '$CONFIG_FILE' not found!"
    exit 1
fi
source "$CONFIG_FILE"
shift   # consume CONFIG_FILE so remaining args are flags

# ---------------------------
# Argument parsing (Step 5 selector)
# ---------------------------
# STEP5_MODE default comes from config; flags override it

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rerun)
      STEP5_MODE="rerun"
      shift
      ;;
    --prod)
      STEP5_MODE="prod"
      shift
      ;;
    --help|-h)
      echo "Usage: $0 CONFIG_FILE [--rerun | --prod]"
      echo
      echo "Step 5 options (override STEP5_MODE in config):"
      echo "  --prod   Run Step 5b"
      echo "  --rerun  Run Step 5a"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      echo "Run with --help to see available options."
      exit 1
      ;;
  esac
done

# ---------------------------
# Config sanity checks
# ---------------------------
: "${LIB_LIST:?LIB_LIST is not set in config}"
: "${WORKDIR:?WORKDIR is not set in config}"
: "${EUK_DIR:?EUK_DIR is not set in config}"
: "${LOG_FILE:?LOG_FILE is not set in config}"
: "${SAMPLE_LIST:?SAMPLE_LIST is not set in config}"
: "${THREADS:?THREADS is not set in config}"
: "${THREADSP:?THREADSP is not set in config}"
: "${QC_DISREGARDED_SUMMARY:?QC_DISREGARDED_SUMMARY is not set in config}"
: "${QC_OUTPUT:?QC_OUTPUT is not set in config}"
: "${QC_STATS_OUTPUT:?QC_STATS_OUTPUT is not set in config}"
: "${SORT_TMPDIR:?SORT_TMPDIR is not set in config}"
: "${STEP5_MODE:?STEP5_MODE is not set in config}"
: "${RSCRIPT_RERUN:?RSCRIPT_RERUN is not set in config}"
: "${RSCRIPT_PROD:?RSCRIPT_PROD is not set in config}"
: "${ALNFILT_MODE:?ALNFILT_MODE is not set in config}"
: "${ALNFILT_MINANI:?ALNFILT_MINANI is not set in config}"
: "${ALNFILT_MAXANI:?ALNFILT_MAXANI is not set in config}"
: "${UNICORN_MINREADS:?UNICORN_MINREADS is not set in config}"
: "${METADMG_SIM_SCORE_LOW:?METADMG_SIM_SCORE_LOW is not set in config}"
: "${METADMG_SIM_SCORE_HIGH:?METADMG_SIM_SCORE_HIGH is not set in config}"
: "${METADMG_HOW_MANY:?METADMG_HOW_MANY is not set in config}"
: "${METADMG_WEIGHT_TYPE:?METADMG_WEIGHT_TYPE is not set in config}"
: "${METADMG_LCA_RANK:?METADMG_LCA_RANK is not set in config}"
: "${METADMG_FIX_NCBI:?METADMG_FIX_NCBI is not set in config}"
: "${METADMG_THREADS:?METADMG_THREADS is not set in config}"
: "${METADMG_SHOWFITS:?METADMG_SHOWFITS is not set in config}"
: "${METADMG_NOPT:?METADMG_NOPT is not set in config}"
: "${METADMG_NBOOTSTRAP:?METADMG_NBOOTSTRAP is not set in config}"
: "${METADMG_DOBOOT:?METADMG_DOBOOT is not set in config}"
: "${METADMG_SEED:?METADMG_SEED is not set in config}"
: "${METADMG_LIB:?METADMG_LIB is not set in config}"

[[ -f "$LIB_LIST" ]] || { echo "[ERROR] LIB_LIST not found: $LIB_LIST" | tee -a "$LOG_FILE"; exit 1; }

: "${MISSING_DATA_LOG:=$WORKDIR/missing_data.log}"

mkdir -p "$(dirname "$LOG_FILE")"
touch "$LOG_FILE"

mkdir -p "$(dirname "$MISSING_DATA_LOG")"
touch "$MISSING_DATA_LOG"

# ---------------------------
# Functions
# ---------------------------
check_success() {
    if [ $? -ne 0 ]; then
        echo "[ERROR] $1 failed. Check the logs for details." | tee -a "$LOG_FILE"
        exit 1
    fi
}

log_step() {
    echo "[$(date)] $1" | tee -a "$LOG_FILE"
}

log_missing() {
    echo "[$(date)] $1" | tee -a "$MISSING_DATA_LOG" >&2
}

# ---------------------------
# Step 1: Build sample.path using find_library.sh
# ---------------------------
log_step "Step 1: Building path list with find_library.sh -> $SAMPLE_LIST"

awk '{ gsub(/\r$/, "", $0); if ($0!="") print $0 }' "$LIB_LIST" \
  | /datasets/caeg_production/_UTILS/find_library.sh \
  | tail -n +2 | awk '$7 != "ARCHIVE"' > "$SAMPLE_LIST"

check_success "Step 1 (build path list)"

log_step "Step 1 done: found $(wc -l < "$SAMPLE_LIST") matching path lines"
if [[ $(wc -l < "$SAMPLE_LIST") != $(wc -l < "$LIB_LIST") ]]; then
  log_step "Warning: Found $(wc -l < "$LIB_LIST") libraries but $(wc -l < "$SAMPLE_LIST") paths."

  awk 'NR==FNR { found[$1]=1; next } FNR>1 { if (!($1 in found)) print $1 }' \
    <(cut -f1 "$SAMPLE_LIST") "$LIB_LIST" \
    | while read -r lib; do
        [[ -n "$lib" ]] && echo "[$(date)] LIBRARY_ID=$lib STEP=1 MESSAGE=No path found in find_library.sh output. Please note, if your data is new the dependencies for ./find_library.sh might not have been updated yet." >> "$MISSING_DATA_LOG"
      done
fi

# ---------------------------
# Step 2: Get meta data from SMDB
# ---------------------------
log_step "Step 2: Fetching metadata from SMDB"

/projects/caeg/apps/fetch-smdb/fetchsmdb.sh "$LIB_LIST" "$WORKDIR/metadata" tsv
check_success "Step 2 (fetch SMDB metadata)"

log_step "Step 2: Validating SMDB metadata output"

metadata_tmp="$(mktemp)"
expected_tmp="$(mktemp)"

awk 'NF > 0 { gsub(/\r$/, "", $1); print $1 }' "$LIB_LIST" | sort -u > "$expected_tmp"

find "$WORKDIR/metadata" -type f -name "*.tsv" -print0 | \
  xargs -0 awk -F'\t' '
    FNR == 1 { next }
    NF > 0 { print }
  ' > "$metadata_tmp" || true

if [[ ! -s "$metadata_tmp" ]]; then
  log_step "Warning: no metadata rows found in $WORKDIR/metadata"
  while read -r lib; do
    [[ -n "$lib" ]] && echo "[$(date)] LIBRARY_ID=$lib STEP=2 MESSAGE=No SMDB metadata rows found at all" >> "$MISSING_DATA_LOG"
  done < "$expected_tmp"
else
  awk -F'\t' '
    BEGIN { OFS="\t" }
    FNR == NR {
      expected[$1] = 1
      next
    }
    {
      lib = $1
      if (lib == "") next
      seen[lib] = 1
      informative = 0
      for (i = 2; i <= NF; i++) {
        val = $i
        gsub(/^[ \t]+|[ \t]+$/, "", val)
        low = tolower(val)
        if (val != "" && low != "na" && low != "nan" && low != "null") {
          informative = 1
          break
        }
      }
      if (informative) {
        has_data[lib] = 1
      }
    }
    END {
      for (lib in expected) {
        if (!(lib in seen)) {
          print "MISSING", lib
        } else if (!(lib in has_data)) {
          print "ALL_NA", lib
        }
      }
    }
  ' "$expected_tmp" "$metadata_tmp" | while IFS=$'\t' read -r status lib; do
    case "$status" in
      MISSING)
        echo "[$(date)] LIBRARY_ID=$lib STEP=2 MESSAGE=Library missing from SMDB metadata output" >> "$MISSING_DATA_LOG"
        ;;
      ALL_NA)
        echo "[$(date)] LIBRARY_ID=$lib STEP=2 MESSAGE=Library present in SMDB metadata output but all metadata fields are NA/empty" >> "$MISSING_DATA_LOG"
        ;;
    esac
  done
fi

rm -f "$metadata_tmp" "$expected_tmp"
log_step "Step 2 done: SMDB metadata validation complete"

# ---------------------------
# Step 3a: Count reads + length stats
# ---------------------------
if [[ -s "$QC_DISREGARDED_SUMMARY" ]] && [[ $(wc -l < "$QC_DISREGARDED_SUMMARY") -gt 1 ]]; then
  log_step "Step 3a skipped: output already exists and is not empty: $QC_DISREGARDED_SUMMARY"
else
  log_step "Step 3a: Counting input total, collapsed total, disregarded and saturated reads + sequence length stats"

  mkdir -p "$(dirname "$QC_DISREGARDED_SUMMARY")"
  mkdir -p "$(dirname "$QC_OUTPUT")"
  mkdir -p "$(dirname "$QC_STATS_OUTPUT")"
  mkdir -p "$SORT_TMPDIR"

  echo -e "library_id\tinput_total_files\tinput_total_reads\tcollapsed_total_files\tcollapsed_total_reads\tdisregarded_files\tdisregarded_reads\tdisregarded_min_len\tdisregarded_avg_len\tdisregarded_max_len\tsaturated_file\tsaturated_reads\tsaturated_min_len\tsaturated_avg_len\tsaturated_max_len" \
    > "$QC_DISREGARDED_SUMMARY"

  parallel -j "$THREADSP" --colsep '\t' '
    lib="{1}"
    base="{4}"

    shopt -s nullglob

    input_total_files=( "$base"/results/reads/raw/Lib_"$lib"_L00*_R1.fastq.gz )
    collapsed_total_files=( "$base"/results/reads/trim/Lib_"$lib"_L00*_collapsed.fastq.gz )
    disregarded_files=( "$base"/results/reads/trim/Lib_"$lib"_L00*_discarded.fastq.gz )
    saturated_file="$base/results/reads/saturated/Lib_${lib}_collapsed.fastq.gz"

    fastq_stats() {
      if (( $# == 0 )); then
        echo -e "NA\tNA\tNA\tNA"
        return 0
      fi
      zcat "$@" 2>/dev/null | awk '"'"'
        NR % 4 == 2 {
          len = length($0)
          n++
          sum += len
          if (min == "" || len < min) min = len
          if (max == "" || len > max) max = len
        }
        END {
          if (n == 0) {
            print "0\tNA\tNA\tNA"
          } else {
            printf "%d\t%d\t%.2f\t%d\n", n, min, sum/n, max
          }
        }
      '"'"'
    }

    fastq_read_count() {
      if (( $# == 0 )); then
        echo "NA"
        return 0
      fi
      zcat "$@" 2>/dev/null | awk '"'"'END { print NR/4 }'"'"'
    }

    input_total_file_list="NA"
    input_total_reads="NA"
    if (( ${#input_total_files[@]} > 0 )); then
      input_total_file_list=$(printf "%s," "${input_total_files[@]}")
      input_total_file_list=${input_total_file_list%,}
      input_total_reads=$(fastq_read_count "${input_total_files[@]}")
    else
      echo "[$(date)] LIBRARY_ID=${lib} STEP=3a MESSAGE=No raw R1 FASTQ files found under $base/results/reads/raw" >> "'"$MISSING_DATA_LOG"'"
    fi

    collapsed_total_file_list="NA"
    collapsed_total_reads="NA"
    if (( ${#collapsed_total_files[@]} > 0 )); then
      collapsed_total_file_list=$(printf "%s," "${collapsed_total_files[@]}")
      collapsed_total_file_list=${collapsed_total_file_list%,}
      collapsed_total_reads=$(fastq_read_count "${collapsed_total_files[@]}")
    else
      echo "[$(date)] LIBRARY_ID=${lib} STEP=3a MESSAGE=No collapsed FASTQ files found under $base/results/reads/trim" >> "'"$MISSING_DATA_LOG"'"
    fi

    disregarded_file_list="NA"
    disregarded_reads="NA"
    disregarded_min="NA"
    disregarded_avg="NA"
    disregarded_max="NA"
    if (( ${#disregarded_files[@]} > 0 )); then
      disregarded_file_list=$(printf "%s," "${disregarded_files[@]}")
      disregarded_file_list=${disregarded_file_list%,}
      read -r disregarded_reads disregarded_min disregarded_avg disregarded_max < <(
        fastq_stats "${disregarded_files[@]}"
      )
    else
      echo "[$(date)] LIBRARY_ID=${lib} STEP=3a MESSAGE=No discarded FASTQ files found under $base/results/reads/trim" >> "'"$MISSING_DATA_LOG"'"
    fi

    saturated_reads="NA"
    saturated_min="NA"
    saturated_avg="NA"
    saturated_max="NA"
    if [[ -e "$saturated_file" ]]; then
      read -r saturated_reads saturated_min saturated_avg saturated_max < <(
        fastq_stats "$saturated_file"
      )
    else
      saturated_file="NA"
      echo "[$(date)] LIBRARY_ID=${lib} STEP=3a MESSAGE=No saturated FASTQ found at $base/results/reads/saturated/Lib_${lib}_collapsed.fastq.gz" >> "'"$MISSING_DATA_LOG"'"
    fi

    echo -e "${lib}\t${input_total_file_list}\t${input_total_reads}\t${collapsed_total_file_list}\t${collapsed_total_reads}\t${disregarded_file_list}\t${disregarded_reads}\t${disregarded_min}\t${disregarded_avg}\t${disregarded_max}\t${saturated_file}\t${saturated_reads}\t${saturated_min}\t${saturated_avg}\t${saturated_max}"
  ' :::: "$SAMPLE_LIST" >> "$QC_DISREGARDED_SUMMARY"

  check_success "Step 3a (QC summary counts)"
  log_step "Step 3a done: QC summary written to $QC_DISREGARDED_SUMMARY"
fi

# ---------------------------
# Step 3b: Get mapped-read counts + seqkit stats
# ---------------------------
log_step "Step 3b: Collecting mapped-read counts + seqkit stats -> $QC_STATS_OUTPUT"

if [[ ! -s "$QC_STATS_OUTPUT" ]] || [[ $(wc -l < "$QC_STATS_OUTPUT") -le 1 ]]; then
  echo -e "library_id\tseqkit_perfect_duplicates\ttotal_no_reads_mapped_uniq\ttotal_no_reads_mapped_uniq_prefilter\t\
R1_num_seqs\tR1_sum_len\tR1_min_len\tR1_avg_len\tR1_max_len\tR1_Q1\tR1_Q2\tR1_Q3\tR1_sum_gap\tR1_N50\tR1_N50_num\tR1_Q20(%)\tR1_Q30(%)\tR1_AvgQual\tR1_GC(%)\tR1_sum_n\t\
R2_num_seqs\tR2_sum_len\tR2_min_len\tR2_avg_len\tR2_max_len\tR2_Q1\tR2_Q2\tR2_Q3\tR2_sum_gap\tR2_N50\tR2_N50_num\tR2_Q20(%)\tR2_Q30(%)\tR2_AvgQual\tR2_GC(%)\tR2_sum_n\t\
R3_num_seqs\tR3_sum_len\tR3_min_len\tR3_avg_len\tR3_max_len\tR3_Q1\tR3_Q2\tR3_Q3\tR3_sum_gap\tR3_N50\tR3_N50_num\tR3_Q20(%)\tR3_Q30(%)\tR3_AvgQual\tR3_GC(%)\tR3_sum_n\t\
R4_num_seqs\tR4_sum_len\tR4_min_len\tR4_avg_len\tR4_max_len\tR4_Q1\tR4_Q2\tR4_Q3\tR4_sum_gap\tR4_N50\tR4_N50_num\tR4_Q20(%)\tR4_Q30(%)\tR4_AvgQual\tR4_GC(%)\tR4_sum_n\t\
microb_in_num_seqs\tmicrob_in_sum_len\tmicrob_in_min_len\tmicrob_in_avg_len\tmicrob_in_max_len\tmicrob_in_Q1\tmicrob_in_Q2\tmicrob_in_Q3\tmicrob_in_sum_gap\tmicrob_in_N50\tmicrob_in_N50_num\tmicrob_in_Q20(%)\tmicrob_in_Q30(%)\tmicrob_in_AvgQual\tmicrob_in_GC(%)\tmicrob_in_sum_n\t\
prefiltered_num_seqs\tprefiltered_sum_len\tprefiltered_min_len\tprefiltered_avg_len\tprefiltered_max_len\tprefiltered_Q1\tprefiltered_Q2\tprefiltered_Q3\tprefiltered_sum_gap\tprefiltered_N50\tprefiltered_N50_num\tprefiltered_Q20(%)\tprefiltered_Q30(%)\tprefiltered_AvgQual\tprefiltered_GC(%)\tprefiltered_sum_n" \
    > "$QC_STATS_OUTPUT"

  process_library() {
    LIB_ID="$1"
    BASE_PATH="$4"

    PERFECT_DUP="$BASE_PATH/logs/reads/derep/Lib_${LIB_ID}_collapsed.log"
    STAT_Z="NA"
    if [[ -f "$PERFECT_DUP" ]]; then
      STAT_Z=$(awk "{print \$2}" "$PERFECT_DUP" | head -n 1)
    else
      echo "[$(date)] LIBRARY_ID=${LIB_ID} STEP=3b MESSAGE=Missing derep log: $PERFECT_DUP" >> "$MISSING_DATA_LOG"
    fi

    MAPPED_READS_UNIQ="$BASE_PATH/results/aligns/merge/Lib_${LIB_ID}_collapsed.bam"
    STAT_G="NA"
    if [[ -f "$MAPPED_READS_UNIQ" ]]; then
      STAT_G=$(samtools view "$MAPPED_READS_UNIQ" | cut -f1 | sort -u --temporary-directory="$SORT_TMPDIR" | wc -l)
    else
      echo "[$(date)] LIBRARY_ID=${LIB_ID} STEP=3b MESSAGE=Missing merged BAM: $MAPPED_READS_UNIQ" >> "$MISSING_DATA_LOG"
    fi

    MAPPED_READS_UNIQ_P="$BASE_PATH/results/prefilter_aligns/merge/Lib_${LIB_ID}_collapsed.bam"
    STAT_H="NA"
    if [[ -f "$MAPPED_READS_UNIQ_P" ]]; then
      STAT_H=$(samtools view "$MAPPED_READS_UNIQ_P" | cut -f1 | sort -u --temporary-directory="$SORT_TMPDIR" | wc -l)
    else
      echo "[$(date)] LIBRARY_ID=${LIB_ID} STEP=3b MESSAGE=Missing prefilter BAM: $MAPPED_READS_UNIQ_P" >> "$MISSING_DATA_LOG"
    fi

    FASTQS=(
      "$BASE_PATH/results/reads/trim/Lib_${LIB_ID}_L001_collapsed.fastq.gz"
      "$BASE_PATH/results/reads/trim/Lib_${LIB_ID}_L002_collapsed.fastq.gz"
      "$BASE_PATH/results/reads/trim/Lib_${LIB_ID}_L003_collapsed.fastq.gz"
      "$BASE_PATH/results/reads/trim/Lib_${LIB_ID}_L004_collapsed.fastq.gz"
      "$BASE_PATH/results/reads/low_complexity/Lib_${LIB_ID}_collapsed.fastq.gz"
      "$BASE_PATH/results/reads/prefilter/extract/Lib_${LIB_ID}_collapsed.fastq.gz"
    )

    SEQKIT_LINE=$(
      parallel -k -j 6 '
        FQ={};
        if [[ -f "$FQ" ]]; then
          seqkit stats -a -T "$FQ" | awk "NR==2 {OFS=\"\t\"; for (i=4; i<=NF; i++) printf \"%s%s\", \$i, (i==NF?ORS:OFS)}"
        else
          echo "[$(date)] LIBRARY_ID='"$LIB_ID"' STEP=3b MESSAGE=Missing seqkit input: $FQ" >> "'"$MISSING_DATA_LOG"'"
          echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
        fi
      ' ::: "${FASTQS[@]}" | paste -sd $'\t' -
    )

    echo -e "${LIB_ID}\t${STAT_Z}\t${STAT_G}\t${STAT_H}\t${SEQKIT_LINE}"
  }

  export -f process_library
  export SORT_TMPDIR
  export MISSING_DATA_LOG

  parallel -j "$THREADSP" --colsep '\t' process_library {1} {2} {3} {4} :::: "$SAMPLE_LIST" >> "$QC_STATS_OUTPUT"
  check_success "Step 3b (mapped reads + seqkit stats)"
  log_step "Step 3b done: QC table with mapped reads + full seqkit stats written to $QC_STATS_OUTPUT"
else
  log_step "Step 3b skipped: $QC_STATS_OUTPUT already exists and is not empty"
fi

# ---------------------------
# Step 5a: Rerun unicorn + metaDMG from production BAMs
# ---------------------------
if [[ "$STEP5_MODE" == "rerun" ]]; then
  log_step "Step 5a selected (--rerun): rerunning unicorn + metaDMG from production BAMs"
  log_step "Step 5a precheck: validating production BAM inputs"

  parallel -j "$THREADSP" --colsep '\t' '
    bam="{4}/results/aligns/merge/Lib_{1}_collapsed.bam"
    if [[ ! -s "$bam" ]]; then
      echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Missing or empty production BAM: $bam" >> "'"$MISSING_DATA_LOG"'"
    fi
  ' :::: "$SAMPLE_LIST"

  log_step "Running unicorn alnfilter..."
  parallel -j "$THREADSP" --colsep '\t' \
      "/projects/caeg/apps/unicorn_300326/unicorn/unicorn alnfilt \
         -b {4}/results/aligns/merge/Lib_{1}_collapsed.bam \
         -t $THREADS --mode $ALNFILT_MODE \
         --outbam $EUK_DIR/{1}.alnfilt.bam \
         --outstat $EUK_DIR/{1}.alnfilt.refstats \
         --minani $ALNFILT_MINANI --maxani $ALNFILT_MAXANI" \
      :::: "$SAMPLE_LIST"
  check_success "Unicorn alnfilt filtering"

  log_step "Running unicorn refstats filter..."
  parallel -j "$THREADSP" --colsep '\t' \
      "/projects/caeg/apps/unicorn_300326/unicorn/unicorn refstats \
         -b {4}/results/aligns/merge/Lib_{1}_collapsed.bam \
         -t $THREADS --minreads $UNICORN_MINREADS --outbam $EUK_DIR/{1}.filtered.bam \
         --outstat $EUK_DIR/{1}.filtered.refstats \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp" \
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
      "/projects/caeg/apps/metaDMG_300326/metaDMG-cpp/metaDMG-cpp lca \
        --names $TAX_PATH_NCBI/taxdump/names.dmp \
        --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
        --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/phylo_norway.acc2taxid.gz) \
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
      "/projects/caeg/apps/metaDMG_300326/metaDMG-cpp/metaDMG-cpp dfit \
        $EUK_DIR/{1}.sort.filtered.bdamage.gz \
        --threads $METADMG_THREADS \
        --names $TAX_PATH_NCBI/taxdump/names.dmp \
        --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
        --showfits $METADMG_SHOWFITS \
        --nopt $METADMG_NOPT \
        --nbootstrap $METADMG_NBOOTSTRAP \
        --seed $METADMG_SEED \
        --lib $METADMG_LIB \
        --out_prefix $EUK_DIR/{1}.sort.filtered" \
      :::: "$SAMPLE_LIST"
  check_success "Damage calculations done"

  log_step "Aggregating lca and dfit metaDMG..."
  parallel -j "$THREADSP" --colsep '\t' \
      "/projects/caeg/apps/metaDMG_300326/metaDMG-cpp/metaDMG-cpp aggregate \
        $EUK_DIR/{1}.sort.filtered.bdamage.gz \
        --names $TAX_PATH_NCBI/taxdump/names.dmp \
        --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
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
      "/projects/caeg/apps/unicorn_300326/unicorn/unicorn taxstats \
         -b $EUK_DIR/{1}.filtered.bam \
         -t $THREADS \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
         --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/phylo_norway.acc2taxid.gz) \
         -k 17 --outstat $EUK_DIR/{1}.filtered.family.taxstats --minreads 3 --rank family" \
      :::: "$SAMPLE_LIST"
    check_success "Unicorn taxstats (family)"

    parallel -j "$THREADSP" --colsep '\t' '
      [[ -s "'"$EUK_DIR"'/{1}.filtered.family.taxstats" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn family taxstats missing or empty" >> "'"$MISSING_DATA_LOG"'"
    ' :::: "$SAMPLE_LIST"

    log_step "Unicorn per taxID statistics (genus level)..."
    parallel -j "$THREADSP" --colsep '\t' \
      "/projects/caeg/apps/unicorn_300326/unicorn/unicorn taxstats \
         -b $EUK_DIR/{1}.filtered.bam \
         -t $THREADS \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
         --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/phylo_norway.acc2taxid.gz) \
         -k 17 --outstat $EUK_DIR/{1}.filtered.genus.taxstats --minreads 3 --rank genus" \
      :::: "$SAMPLE_LIST"
    check_success "Unicorn taxstats (genus)"

    parallel -j "$THREADSP" --colsep '\t' '
      [[ -s "'"$EUK_DIR"'/{1}.filtered.genus.taxstats" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn genus taxstats missing or empty" >> "'"$MISSING_DATA_LOG"'"
    ' :::: "$SAMPLE_LIST"

    log_step "Unicorn per taxID statistics (species level)..."
    parallel -j "$THREADSP" --colsep '\t' \
      "/projects/caeg/apps/unicorn_300326/unicorn/unicorn taxstats \
         -b $EUK_DIR/{1}.filtered.bam \
         -t $THREADS \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
         --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/phylo_norway.acc2taxid.gz) \
         -k 17 --outstat $EUK_DIR/{1}.filtered.spec.taxstats --minreads 3 --rank species" \
      :::: "$SAMPLE_LIST"
    check_success "Unicorn taxstats (species)"

    parallel -j "$THREADSP" --colsep '\t' '
      [[ -s "'"$EUK_DIR"'/{1}.filtered.spec.taxstats" ]] || echo "[$(date)] LIBRARY_ID={1} STEP=5a MESSAGE=Unicorn species taxstats missing or empty" >> "'"$MISSING_DATA_LOG"'"
    ' :::: "$SAMPLE_LIST"


  log_step "Running Rscript for merging and extracting TSV tables on species and genus level"
  Rscript "$RSCRIPT_RERUN" "$CONFIG_FILE"
  check_success "Rscript complete."
fi

# ---------------------------
# Step 5b: Use production unicorn/metaDMG outputs
# ---------------------------
if [[ "$STEP5_MODE" == "prod" ]]; then
  log_step "Step 5b selected (--prod): using production unicorn/metaDMG outputs"

  log_step "Merging TSVs per library (keep one header) and filtering col4 > $UNICORN_MINREADS..."
  parallel -j "$THREADSP" --colsep '\t' '
    in_dir="{4}/stats/shards/unicorn/";
    out="'"$EUK_DIR"'/{1}.filtered.unicorn.refstats.tsv";
    mkdir -p "$(dirname "$out")"
    shopt -s nullglob
    files=("$in_dir"/*.tsv)

    if (( ${#files[@]} == 0 )); then
      echo "[WARN] No TSVs found for {1} in $in_dir" >&2
      echo "[$(date)] LIBRARY_ID={1} STEP=5b MESSAGE=No unicorn shard TSVs found in $in_dir" >> "'"$MISSING_DATA_LOG"'"
      exit 0
    fi

    awk -F"\t" '"'"'FNR==1 { if (NR==1) print; next } $4+0 > '"$UNICORN_MINREADS"' { print }'"'"' "${files[@]}" > "$out"
    if [[ ! -s "$out" ]]; then
      echo "[$(date)] LIBRARY_ID={1} STEP=5b MESSAGE=Merged unicorn refstats output missing or empty: $out" >> "'"$MISSING_DATA_LOG"'"
    fi
  ' :::: "$SAMPLE_LIST"
  check_success "TSV merge+filter"

  log_step "Symlinking metaDMG aggregate stat files into EUK_DIR..."
  parallel -j "$THREADSP" --colsep '\t' \
    'src="{4}/results/metadmg/aggregate/Lib_{1}_collapsed.stat.gz";
     dst="'"$EUK_DIR"'/{1}.collapsed.stat.gz";

     if [[ ! -e "$src" ]]; then
       echo "[WARN] Missing: $src" >&2
       echo "[$(date)] LIBRARY_ID={1} STEP=5b MESSAGE=Missing production metaDMG stat.gz: $src" >> "'"$MISSING_DATA_LOG"'"
       exit 0
     fi

     ln -sfn "$src" "$dst"' \
    :::: "$SAMPLE_LIST"
  check_success "Symlink metaDMG aggregate stats"

  log_step "Running Rscript for merging and extracting TSV tables on species and genus level"
  Rscript "$RSCRIPT_PROD" "$CONFIG_FILE"
  check_success "Rscript complete."
fi

log_step "Butteracid is done."