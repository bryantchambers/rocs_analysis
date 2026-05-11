#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"

# Defaults are "final-like" and intentionally heavy.
N_BOOT="${N_BOOT:-120}"
N_PERM="${N_PERM:-700}"
TOP_N="${TOP_N:-5}"

"${RSCRIPT}" "${S}/06_full_eval_top5.R" "--n_boot=${N_BOOT}" "--n_perm=${N_PERM}" "--top_n=${TOP_N}"
"${RSCRIPT}" "${S}/07_full_eval_report.R"

echo "[networkQC] full eval complete"

