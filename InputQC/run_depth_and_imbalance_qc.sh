#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

echo "[InputQC] Running rarefaction/depth QC"
Rscript InputQC/rarefaction_depth_qc/scripts/01_rarefaction_depth_qc.R

echo "[InputQC] Running core/age imbalance QC"
Rscript InputQC/core_age_imbalance/scripts/01_core_age_imbalance.R

echo "[InputQC] Complete"
