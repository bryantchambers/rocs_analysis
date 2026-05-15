#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"

# This is a focused sensitivity run, not a final preservation run.
N_PERM="${N_PERM:-100}"

"${RSCRIPT}" "${S}/01_data_prep_input_variants.R"
"${RSCRIPT}" "${S}/02_eval_input_sensitivity.R" "--n_perm=${N_PERM}"
"${RSCRIPT}" "${S}/04_input_decision_review.R"

echo "[InputQC] complete"
