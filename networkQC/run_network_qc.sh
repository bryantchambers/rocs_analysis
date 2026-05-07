#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"

"${RSCRIPT}" "${S}/00_config_check.R"
"${RSCRIPT}" "${S}/01_wgcna_qc_basics.R"
"${RSCRIPT}" "${S}/02_wgcna_parameter_sweep.R"
"${RSCRIPT}" "${S}/03_leiden_modules.R"
"${RSCRIPT}" "${S}/04_compare_wgcna_leiden.R"

echo "[networkQC] complete"

