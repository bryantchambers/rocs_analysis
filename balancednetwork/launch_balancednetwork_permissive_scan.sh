#!/usr/bin/env bash
set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${BASE}/results/logs"
mkdir -p "${LOG_DIR}"

RUN_ID="${BALANCED_RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
OUT_LOG="${LOG_DIR}/permissive_scan_${RUN_ID}.out"
PID_FILE="${LOG_DIR}/permissive_scan_${RUN_ID}.pid"
PROGRESS_TSV="${BALANCED_RUN_PROGRESS_TSV:-${LOG_DIR}/permissive_scan_${RUN_ID}.progress.tsv}"

export BALANCED_RUN_ID="${RUN_ID}"
export BALANCED_RUN_PROGRESS_TSV="${PROGRESS_TSV}"

setsid nohup bash "${BASE}/run_balancednetwork_permissive_scan.sh" "$@" >"${OUT_LOG}" 2>&1 < /dev/null &
PID=$!
printf '%s\n' "${PID}" >"${PID_FILE}"

printf 'run_id=%s\npid=%s\npid_file=%s\nprogress_tsv=%s\nout_log=%s\n' \
  "${RUN_ID}" \
  "${PID}" \
  "${PID_FILE}" \
  "${PROGRESS_TSV}" \
  "${OUT_LOG}"
