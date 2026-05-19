#!/usr/bin/env bash
# Download the prebuilt SnpEff database for Sorghum bicolor.
#
# Uses `snpEff download` to fetch the prebuilt database from the SnpEff
# repository rather than building from local FASTA/GFF files.
#
# Usage (from any directory within the repo):
#   bash analysis/scripts/snpeff_prep.sh
#
# To verify the exact genome name available in the SnpEff repo:
#   bash docker/run.sh snpEff databases | grep -i sorghum
set -euo pipefail

THESIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SNPEFF_DIR="${THESIS_DIR}/resources/snpeff"
LOG_DIR="${THESIS_DIR}/analysis/logs"
mkdir -p "${LOG_DIR}"
LOG="${LOG_DIR}/snpeff_prep.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

DB="Sorghum_bicolor"

DATA_DIR="${SNPEFF_DIR}/data"
CONFIG="${SNPEFF_DIR}/snpEff.config"

mkdir -p "${DATA_DIR}"

# Minimal config so snpEff knows where to place the downloaded database.
cat > "${CONFIG}" <<EOF
data.dir = ${SNPEFF_DIR}/data
EOF

echo "[snpeff_prep] downloading prebuilt SnpEff database for ${DB}"
snpEff download \
    -dataDir "${DATA_DIR}" \
    -v \
    "${DB}"

echo "[snpeff_prep] done — database written to ${DATA_DIR}/${DB}"
