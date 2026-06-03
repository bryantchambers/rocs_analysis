#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"

"${RSCRIPT}" "${S}/04_balanced_parameter_sweep.R"
"${RSCRIPT}" "${S}/05_balanced_decision_matrix.R"

echo "[balancednetwork] qc sweep+decision complete"
