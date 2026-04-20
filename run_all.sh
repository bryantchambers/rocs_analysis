#!/usr/bin/env bash
# run_all.sh — ROCS WGCNA + HMM + EMP/TEA pipeline
#
# Usage:
#   bash run_all.sh [--start STEP]
#
# Steps:
#   1  Data prep (filter, VST)
#   2  WGCNA consensus modules
#   3  HMM ecological states
#   4  EMP/SAP metabolic potential
#   5  TEA indices vs EMP comparison

set -euo pipefail

RSCRIPT="/maps/projects/fernandezguerra/apps/opt/conda/envs/r/bin/Rscript"
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${PROJ}/scripts"

START=1
while [[ $# -gt 0 ]]; do
  case $1 in
    --start) START="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

log() { echo "[$(date '+%H:%M:%S')] $*"; }

(( START <= 1 )) && { log "=== 1. Data prep ===";    "${RSCRIPT}" "${S}/01_data_prep.R"; }
(( START <= 2 )) && { log "=== 2. WGCNA ===";        "${RSCRIPT}" "${S}/02_wgcna.R"; }
(( START <= 3 )) && { log "=== 3. HMM states ===";   "${RSCRIPT}" "${S}/03_hmm_states.R"; }
(( START <= 4 )) && { log "=== 4. EMP/SAP ===";      "${RSCRIPT}" "${S}/04_emp.R"; }
(( START <= 5 )) && { log "=== 5. TEA vs EMP ===";   "${RSCRIPT}" "${S}/05_tea_vs_emp.R"; }

log "=== Done ==="
