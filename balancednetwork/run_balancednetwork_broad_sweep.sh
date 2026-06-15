#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"
LOG_DIR="${BASE}/results/logs"
mkdir -p "${LOG_DIR}"

if [[ -z "${BALANCED_BROAD_SCAN_RUN_ID:-}" ]]; then
  export BALANCED_BROAD_SCAN_RUN_ID="$(date +%Y%m%d_%H%M%S)"
fi

if [[ -z "${BALANCED_BROAD_SWEEP_PROGRESS_TSV:-}" ]]; then
  export BALANCED_BROAD_SWEEP_PROGRESS_TSV="${LOG_DIR}/broad_scan_${BALANCED_BROAD_SCAN_RUN_ID}.progress.tsv"
fi

export BALANCED_BROAD_RECUT_TIMEOUT_SEC="${BALANCED_BROAD_RECUT_TIMEOUT_SEC:-600}"
export BALANCED_WGCNA_THREADS="${BALANCED_WGCNA_THREADS:-8}"

echo "[balancednetwork] broad scan run_id=${BALANCED_BROAD_SCAN_RUN_ID}"
echo "[balancednetwork] progress=${BALANCED_BROAD_SWEEP_PROGRESS_TSV}"
echo "[balancednetwork] recut_timeout_sec=${BALANCED_BROAD_RECUT_TIMEOUT_SEC}"
echo "[balancednetwork] wgcna_threads=${BALANCED_WGCNA_THREADS}"

"${RSCRIPT}" "${S}/04b_balanced_broad_power_scan.R"
"${RSCRIPT}" "${S}/04c_balanced_broad_recut_sweep.R"
"${RSCRIPT}" "${S}/04d_balanced_broad_candidates.R"

echo "[balancednetwork] broad sweep complete"
