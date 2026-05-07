#!/usr/bin/env bash
# run_all.sh — ROCS WGCNA, HMM, EMP/TEA, driver, and story pipeline
#
# Usage:
#   bash run_all.sh [--start STEP] [--mode build|final]
#
# Steps:
#   01   Data prep (filter, CLR/VST object)
#   02   WGCNA consensus modules
  #   02b  WGCNA stability diagnostics
#   03   HMM ecological states
#   04   EMP/SAP metabolic potential
#   05   TEA indices vs EMP comparison
#   06   Diagnostic comparison figures
#   06b  Bryant figure set
#   07   Taxon importance (kME + RF)
#   07b  Taxon importance (FuzzyForest)
#   08   Network statistics
#   09   Driver integration
#   10   Climate sensitivity
#   11   State-specific networks
#   12   Functional linkage
#   13   State transition meta-network
#   14   Driver quadrants
#   15   State functional breakdown
#   16   Final story visualization

set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${PROJ}/scripts"

START="01"
MODE=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --start) START="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

log() { echo "[$(date '+%H:%M:%S')] $*"; }

normalize_step() {
  local step="$1"
  if [[ "${step}" =~ ^0*([0-9]+)([[:alpha:]]*)$ ]]; then
    printf "%02d%s" "$((10#${BASH_REMATCH[1]}))" "${BASH_REMATCH[2]}"
  else
    printf "%s" "${step}"
  fi
}

START="$(normalize_step "${START}")"

PIPELINE=(
  "01|Data prep|01_data_prep.R"
  "02|WGCNA consensus modules|02_wgcna.R"
  "02b|WGCNA stability diagnostics|02b_wgcna_stability.R"
  "03|HMM ecological states|03_hmm_states.R"
  "04|EMP/SAP metabolic potential|04_emp.R"
  "05|TEA indices vs EMP comparison|05_tea_vs_emp.R"
  "06|Diagnostic comparison figures|06_figures.R"
  "06b|Bryant figure set|06b_bryantfigures.R"
  "07|Taxon importance (kME + RF)|07_taxon_importance.R"
  "07b|Taxon importance (FuzzyForest)|07b_taxon_importance_fuzzy.R"
  "08|Network statistics|08_network_statistics.R"
  "09|Driver integration|09_driver_integration.R"
  "10|Climate sensitivity|10_climate_sensitivity.R"
  "11|State-specific networks|11_state_networks.R"
  "12|Functional linkage|12_functional_linkage.R"
  "13|State transition meta-network|13_state_transition_network.R"
  "14|Driver quadrants|14_driver_quadrants.R"
  "15|State functional breakdown|15_state_functional_breakdown.R"
  "16|Final story visualization|16_final_story_visualization.R"
)

start_index=-1
for i in "${!PIPELINE[@]}"; do
  IFS="|" read -r step_id _ _ <<< "${PIPELINE[$i]}"
  if [[ "${step_id}" == "${START}" ]]; then
    start_index="${i}"
    break
  fi
done

if [[ "${start_index}" -lt 0 ]]; then
  echo "Unknown --start step: ${START}" >&2
  echo "Known steps:" >&2
  for entry in "${PIPELINE[@]}"; do
    IFS="|" read -r step_id step_name _ <<< "${entry}"
    echo "  ${step_id}  ${step_name}" >&2
  done
  exit 1
fi

if [[ "${RSCRIPT}" == */* ]]; then
  [[ -x "${RSCRIPT}" ]] || { echo "Rscript is not executable: ${RSCRIPT}" >&2; exit 1; }
else
  command -v "${RSCRIPT}" >/dev/null 2>&1 || { echo "Rscript not found on PATH: ${RSCRIPT}" >&2; exit 1; }
fi

for i in "${!PIPELINE[@]}"; do
  (( i < start_index )) && continue
  IFS="|" read -r step_id step_name script_name <<< "${PIPELINE[$i]}"
  script_path="${S}/${script_name}"
  [[ -f "${script_path}" ]] || { echo "Pipeline script not found: ${script_path}" >&2; exit 1; }
  log "=== ${step_id}. ${step_name} ==="
  if [[ -n "${MODE}" && ( "${step_id}" == "02" || "${step_id}" == "02b" ) ]]; then
    "${RSCRIPT}" "${script_path}" "--mode=${MODE}"
  else
    "${RSCRIPT}" "${script_path}"
  fi
done

log "=== Done ==="
