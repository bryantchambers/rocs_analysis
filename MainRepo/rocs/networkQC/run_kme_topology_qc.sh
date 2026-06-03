#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${BASE}/results/kme_topology_run.log"

mkdir -p "${BASE}/results"

{
  echo "[networkQC] kME/topology QC started: $(date)"
  "${RSCRIPT}" "${BASE}/scripts/10_kme_topology_review.R"
  echo "[networkQC] kME/topology QC complete: $(date)"
} 2>&1 | tee "${LOG}"
