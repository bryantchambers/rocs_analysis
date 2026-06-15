#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"
LOG_DIR="${BASE}/results/logs"
mkdir -p "${LOG_DIR}"

if [[ -z "${BALANCED_NEIGHBORHOOD_RUN_ID:-}" ]]; then
  export BALANCED_NEIGHBORHOOD_RUN_ID="$(date +%Y%m%d_%H%M%S)"
fi

if [[ -z "${BALANCED_NEIGHBORHOOD_PROGRESS_TSV:-}" ]]; then
  export BALANCED_NEIGHBORHOOD_PROGRESS_TSV="${LOG_DIR}/neighborhood_scan_${BALANCED_NEIGHBORHOOD_RUN_ID}.progress.tsv"
fi

export BALANCED_WGCNA_THREADS="${BALANCED_WGCNA_THREADS:-8}"

echo "[balancednetwork] neighborhood scan run_id=${BALANCED_NEIGHBORHOOD_RUN_ID}"
echo "[balancednetwork] progress=${BALANCED_NEIGHBORHOOD_PROGRESS_TSV}"
echo "[balancednetwork] wgcna_threads=${BALANCED_WGCNA_THREADS}"

"${RSCRIPT}" "${S}/04e_balanced_neighborhood_scan.R"
"${RSCRIPT}" "${S}/04f_balanced_neighborhood_candidates.R"

echo "[balancednetwork] neighborhood scan complete"
