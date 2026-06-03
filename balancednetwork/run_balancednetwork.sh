#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${BASE}/scripts"

MODE="build"
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
  N_BOOT_DEFAULT="1000"
  N_PERM_DEFAULT="700"
else
  N_BOOT_DEFAULT="100"
  N_PERM_DEFAULT="200"
fi

N_BOOT="${N_BOOT:-$N_BOOT_DEFAULT}"
N_PERM="${N_PERM:-$N_PERM_DEFAULT}"
BIN_WIDTH_KYR="${BIN_WIDTH_KYR:-10}"
FORCE="${FORCE:-0}"

"${RSCRIPT}" "${S}/00_balance_design.R" "--bin_width_kyr=${BIN_WIDTH_KYR}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/01_balanced_wgcna_exp3.R" "--mode=${MODE}" "--n_perm=${N_PERM}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/02_balanced_wgcna_stability.R" "--mode=${MODE}" "--n_boot=${N_BOOT}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/03_balanced_network_qc_report.R" "--mode=${MODE}"

echo "[balancednetwork] complete: mode=${MODE}, n_boot=${N_BOOT}, n_perm=${N_PERM}, bin_width=${BIN_WIDTH_KYR}"
