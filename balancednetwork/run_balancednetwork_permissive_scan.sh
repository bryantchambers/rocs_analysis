#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"
LOG_DIR="${BASE}/results/logs"
mkdir -p "${LOG_DIR}"

MODE="build"
for a in "$@"; do
  case "$a" in
    --mode=*) MODE="${a#--mode=}" ;;
  esac
done

if [[ "$MODE" != "build" && "$MODE" != "final" ]]; then
  echo "invalid mode: $MODE (use --mode=build or --mode=final)" >&2
  exit 1
fi

if [[ "$MODE" == "final" ]]; then
  N_BOOT_DEFAULT="1000"
  N_PERM_DEFAULT="700"
else
  N_BOOT_DEFAULT="100"
  N_PERM_DEFAULT="200"
fi

RUN_ID="${BALANCED_RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
STAGE1_VARIANT="${BALANCED_STAGE1_VARIANT:-permissive_abund10}"
N_BOOT="${N_BOOT:-$N_BOOT_DEFAULT}"
N_PERM="${N_PERM:-$N_PERM_DEFAULT}"
BIN_WIDTH_KYR="${BIN_WIDTH_KYR:-10}"
FORCE="${FORCE:-0}"
BALANCED_NEIGHBORHOOD_EVAL_TOP_N="${BALANCED_NEIGHBORHOOD_EVAL_TOP_N:-4}"

export BALANCED_RUN_ID="${RUN_ID}"
export BALANCED_STAGE1_VARIANT="${STAGE1_VARIANT}"
export BALANCED_BROAD_SCAN_RUN_ID="${BALANCED_BROAD_SCAN_RUN_ID:-${RUN_ID}}"
export BALANCED_NEIGHBORHOOD_RUN_ID="${BALANCED_NEIGHBORHOOD_RUN_ID:-${RUN_ID}}"
export BALANCED_RUN_PROGRESS_TSV="${BALANCED_RUN_PROGRESS_TSV:-${LOG_DIR}/permissive_scan_${RUN_ID}.progress.tsv}"
export BALANCED_BROAD_SWEEP_PROGRESS_TSV="${BALANCED_BROAD_SWEEP_PROGRESS_TSV:-${BALANCED_RUN_PROGRESS_TSV}}"
export BALANCED_NEIGHBORHOOD_PROGRESS_TSV="${BALANCED_NEIGHBORHOOD_PROGRESS_TSV:-${BALANCED_RUN_PROGRESS_TSV}}"
export BALANCED_SWEEP_PROGRESS_TSV="${BALANCED_SWEEP_PROGRESS_TSV:-${BALANCED_RUN_PROGRESS_TSV}}"
export BAL_PROGRESS_TSV="${BAL_PROGRESS_TSV:-${BALANCED_RUN_PROGRESS_TSV}}"
export BALANCED_NEIGHBORHOOD_EVAL_TOP_N="${BALANCED_NEIGHBORHOOD_EVAL_TOP_N}"
export BALANCED_BROAD_RECUT_TIMEOUT_SEC="${BALANCED_BROAD_RECUT_TIMEOUT_SEC:-600}"
export BALANCED_WGCNA_THREADS="${BALANCED_WGCNA_THREADS:-8}"

echo "[balancednetwork] permissive scan run_id=${BALANCED_RUN_ID}"
echo "[balancednetwork] stage1_variant=${BALANCED_STAGE1_VARIANT}"
echo "[balancednetwork] progress=${BALANCED_RUN_PROGRESS_TSV}"
echo "[balancednetwork] neighborhood_top_n=${BALANCED_NEIGHBORHOOD_EVAL_TOP_N}"
echo "[balancednetwork] wgcna_threads=${BALANCED_WGCNA_THREADS}"

"${RSCRIPT}" "${S}/00_balanced_input_prep.R"
"${RSCRIPT}" "${S}/00_balance_design.R" "--bin_width_kyr=${BIN_WIDTH_KYR}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/01_balanced_wgcna_exp3.R" "--mode=${MODE}" "--n_perm=${N_PERM}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/02_balanced_wgcna_stability.R" "--mode=${MODE}" "--n_boot=${N_BOOT}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/03_balanced_network_qc_report.R" "--mode=${MODE}"
"${RSCRIPT}" "${S}/04b_balanced_broad_power_scan.R"
"${RSCRIPT}" "${S}/04c_balanced_broad_recut_sweep.R"
"${RSCRIPT}" "${S}/04d_balanced_broad_candidates.R"
"${RSCRIPT}" "${S}/04e_balanced_neighborhood_scan.R"
"${RSCRIPT}" "${S}/04f_balanced_neighborhood_candidates.R"

echo "[balancednetwork] permissive scan complete"
