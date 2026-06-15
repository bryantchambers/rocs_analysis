#!/usr/bin/env bash
set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${BASE}/results/logs"
mkdir -p "${LOG_DIR}"

RUN_ID="${BALANCED_BROAD_SCAN_RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
OUT_LOG="${LOG_DIR}/broad_scan_${RUN_ID}.out"
PID_FILE="${LOG_DIR}/broad_scan_${RUN_ID}.pid"
PROGRESS_TSV="${LOG_DIR}/broad_scan_${RUN_ID}.progress.tsv"

export BALANCED_BROAD_SCAN_RUN_ID="${RUN_ID}"
export BALANCED_BROAD_SWEEP_PROGRESS_TSV="${BALANCED_BROAD_SWEEP_PROGRESS_TSV:-${PROGRESS_TSV}}"
export BALANCED_BROAD_RECUT_TIMEOUT_SEC="${BALANCED_BROAD_RECUT_TIMEOUT_SEC:-600}"
export BALANCED_WGCNA_THREADS="${BALANCED_WGCNA_THREADS:-8}"

setsid nohup bash "${BASE}/run_balancednetwork_broad_sweep.sh" >"${OUT_LOG}" 2>&1 < /dev/null &
echo $! > "${PID_FILE}"

printf 'run_id=%s\npid=%s\npid_file=%s\nprogress_tsv=%s\nout_log=%s\n' \
  "${RUN_ID}" \
  "$(cat "${PID_FILE}")" \
  "${PID_FILE}" \
  "${BALANCED_BROAD_SWEEP_PROGRESS_TSV}" \
  "${OUT_LOG}"
