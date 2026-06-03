#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"

MODE="final"
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
  N_BOOT_DEFAULT="500"
  N_PERM_DEFAULT="700"
else
  N_BOOT_DEFAULT="100"
  N_PERM_DEFAULT="200"
fi

N_BOOT="${N_BOOT:-$N_BOOT_DEFAULT}"
N_PERM="${N_PERM:-$N_PERM_DEFAULT}"
TOP_N="${TOP_N:-5}"
FORCE="${FORCE:-0}"

"${RSCRIPT}" "${S}/06_balanced_full_eval_top5.R" "--mode=${MODE}" "--n_boot=${N_BOOT}" "--n_perm=${N_PERM}" "--top_n=${TOP_N}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/07_balanced_full_eval_report.R"

echo "[balancednetwork] full eval complete: mode=${MODE}, n_boot=${N_BOOT}, n_perm=${N_PERM}, top_n=${TOP_N}"
