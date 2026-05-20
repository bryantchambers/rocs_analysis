# ---------------------------
# Paths
# ---------------------------
LIB_LIST="/projects/caeg/scratch/nare_data/YW_data/sample.list.lib"
PROD_PATH="/datasets/caeg_production/libraries/"
EUK_DIR="/projects/caeg/scratch/nare_data/YW_data/eukaryotic"
WORKDIR="/projects/caeg/scratch/nare_data/YW_data/"
LOG_FILE="/projects/caeg/scratch/nare_data/YW_data/pipeline.log"
SAMPLE_LIST="/projects/caeg/scratch/nare_data/YW_data/sample.path"
LOG_DIR="/projects/caeg/scratch/nare_data/YW_data/logs"
DB_PATH="/datasets/caeg_dataset/references/ncbi/20250530/data/wgs_eukaryota"
DB_PATH_clean="/datasets/caeg_dataset/references/ncbi/20250530/data/"
TAX_PATH_NCBI="/datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/"

# ---------------------------
# Threading
# ---------------------------
THREADS="10"
THREADSP="5"

# ---------------------------
# Step 3: QC output paths
# ---------------------------
QC_DISREGARDED_SUMMARY="${WORKDIR}/disregarded_reads_counts.tsv"
QC_OUTPUT="${WORKDIR}/qc_metadata.filtered.tsv"
QC_STATS_OUTPUT="${WORKDIR}/QC_stats.tsv"
SORT_TMPDIR="/projects/caeg/scratch/nare_data/YW_data/tmp"

# ---------------------------
# Step 5: Mode default (overridden by --rerun / --prod flags)
# ---------------------------
STEP5_MODE="prod"

# ---------------------------
# Step 5: R scripts
# ---------------------------
RSCRIPT_RERUN="butteracid_rerun.R"
RSCRIPT_PROD="butteracid_production.R"

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