#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"
LOG_DIR="${BASE}/results/logs"
mkdir -p "${LOG_DIR}"

if [[ -z "${BALANCED_BROAD_SWEEP_PROGRESS_TSV:-}" ]]; then
  export BALANCED_BROAD_SWEEP_PROGRESS_TSV="${LOG_DIR}/broad_sweep_progress_$(date +%Y%m%d_%H%M%S).tsv"
fi

"${RSCRIPT}" "${S}/04b_balanced_broad_power_scan.R"
"${RSCRIPT}" "${S}/04c_balanced_broad_recut_sweep.R"
"${RSCRIPT}" "${S}/04d_balanced_broad_candidates.R"

echo "[balancednetwork] broad sweep complete"
