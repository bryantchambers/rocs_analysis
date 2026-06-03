#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

export WGCNA_HMM_RESULTS_SUFFIX="${WGCNA_HMM_RESULTS_SUFFIX:-balanced_default_final_$(date +%Y%m%d)}"
export WGCNA_HMM_INPUT_STRATEGY="${WGCNA_HMM_INPUT_STRATEGY:-balanced}"
export WGCNA_HMM_WGCNA_PROFILE="${WGCNA_HMM_WGCNA_PROFILE:-balanced_top3}"
export WGCNA_HMM_INPUT_TAX_DAMAGE="${WGCNA_HMM_INPUT_TAX_DAMAGE:-/src/results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz}"
export WGCNA_HMM_INPUT_KEGG_MODS="${WGCNA_HMM_INPUT_KEGG_MODS:-/src/data/functional/kegg-modules-summary-rocs.tsv.gz}"

RESULTS_DIR="results/wgcna_hmm_${WGCNA_HMM_RESULTS_SUFFIX}"
LOG_DIR="${RESULTS_DIR}/logs"
mkdir -p "${LOG_DIR}"

{
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Bryant merge final workflow"
  echo "repo=$(pwd)"
  echo "results=${RESULTS_DIR}"
  echo "input_strategy=${WGCNA_HMM_INPUT_STRATEGY}"
  echo "wgcna_profile=${WGCNA_HMM_WGCNA_PROFILE}"
  echo "hmm_input_balancing_status=pending_review_with_M"
  echo "tax_damage=${WGCNA_HMM_INPUT_TAX_DAMAGE}"
  echo "kegg_mods=${WGCNA_HMM_INPUT_KEGG_MODS}"
  echo "command=Rscript code/wgcna_hmm/run_workflow.R --mode=final"
} >> "${LOG_DIR}/overnight_run.log"

Rscript code/wgcna_hmm/run_workflow.R --mode=final >> "${LOG_DIR}/overnight_run.log" 2>&1
