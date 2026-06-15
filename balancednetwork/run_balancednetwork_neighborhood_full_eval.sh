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
  N_BOOT_DEFAULT="500"
  N_PERM_DEFAULT="700"
else
  N_BOOT_DEFAULT="100"
  N_PERM_DEFAULT="200"
fi

N_BOOT="${N_BOOT:-$N_BOOT_DEFAULT}"
N_PERM="${N_PERM:-$N_PERM_DEFAULT}"
TOP_N="${TOP_N:-8}"
FORCE="${FORCE:-0}"

if [[ -z "${BALANCED_FULL_EVAL_SOURCE_TSV:-}" ]]; then
  if [[ -n "${BALANCED_NEIGHBORHOOD_RUN_ID:-}" ]]; then
    export BALANCED_FULL_EVAL_SOURCE_TSV="${BASE}/results/qc/neighborhood_scan/${BALANCED_NEIGHBORHOOD_RUN_ID}/neighborhood_settings_to_full_eval.tsv"
    export BALANCED_FULL_EVAL_DIR="${BASE}/results/qc/neighborhood_scan/${BALANCED_NEIGHBORHOOD_RUN_ID}/full_eval"
    export BALANCED_FULL_EVAL_REPORT="${BASE}/results/qc/neighborhood_scan/${BALANCED_NEIGHBORHOOD_RUN_ID}/BALANCED_NEIGHBORHOOD_FULL_EVAL_REPORT.md"
    export BALANCED_FULL_EVAL_RANKED_TSV="${BASE}/results/qc/neighborhood_scan/${BALANCED_NEIGHBORHOOD_RUN_ID}/full_eval/all_settings_ranked.tsv"
  else
    export BALANCED_FULL_EVAL_SOURCE_TSV="${BASE}/results/qc/neighborhood_scan/neighborhood_settings_to_full_eval.tsv"
    export BALANCED_FULL_EVAL_DIR="${BASE}/results/qc/neighborhood_scan/full_eval"
    export BALANCED_FULL_EVAL_REPORT="${BASE}/results/qc/neighborhood_scan/BALANCED_NEIGHBORHOOD_FULL_EVAL_REPORT.md"
    export BALANCED_FULL_EVAL_RANKED_TSV="${BASE}/results/qc/neighborhood_scan/full_eval/all_settings_ranked.tsv"
  fi
fi

if [[ -z "${BAL_PROGRESS_TSV:-}" ]]; then
  export BAL_PROGRESS_TSV="${BASE}/results/logs/neighborhood_full_eval_${BALANCED_NEIGHBORHOOD_RUN_ID:-default}.progress.tsv"
fi

echo "[balancednetwork] full-eval progress=${BAL_PROGRESS_TSV}"

"${RSCRIPT}" "${S}/06_balanced_full_eval_top5.R" "--mode=${MODE}" "--n_boot=${N_BOOT}" "--n_perm=${N_PERM}" "--top_n=${TOP_N}" "--force=${FORCE}"
"${RSCRIPT}" "${S}/07_balanced_full_eval_report.R"

echo "[balancednetwork] neighborhood full eval complete: mode=${MODE}, n_boot=${N_BOOT}, n_perm=${N_PERM}, top_n=${TOP_N}"
