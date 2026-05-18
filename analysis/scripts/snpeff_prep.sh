#!/usr/bin/env bash
# Build the SnpEff database for Sorghum bicolor NCBIv3.
#
# Copies the reference genome FASTA and GFF annotation into the SnpEff data
# directory, writes a snpEff.config, then runs `snpEff build` via Docker.
#
# Usage (from any directory within the repo):
#   bash analysis/scripts/snpeff_prep.sh
set -euo pipefail

THESIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SNPEFF_DIR="${THESIS_DIR}/resources/snpeff"
RUN="${THESIS_DIR}/docker/run.sh"

GENOME="Sorghum_bicolor_NCBIv3"
REF_FA="${THESIS_DIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
ANNOT_GFF="${THESIS_DIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff"

DATA_DIR="${SNPEFF_DIR}/data/${GENOME}"
CONFIG="${SNPEFF_DIR}/snpEff.config"

mkdir -p "${DATA_DIR}"

# --- snpEff.config ---------------------------------------------------------
cat > "${CONFIG}" <<EOF
data.dir = ${SNPEFF_DIR}/data
${GENOME}.genome : Sorghum bicolor NCBIv3
${GENOME}
EOF

# --- genome inputs ---------------------------------------------------------
echo "[snpeff_prep] copying genome FASTA → ${DATA_DIR}/sequences.fa"
cp "${REF_FA}" "${DATA_DIR}/sequences.fa"

echo "[snpeff_prep] copying annotation GFF → ${DATA_DIR}/genes.gff"
cp "${ANNOT_GFF}" "${DATA_DIR}/genes.gff"

# --- build -----------------------------------------------------------------
echo "[snpeff_prep] running snpEff build for ${GENOME}"
bash "${RUN}" snpEff build \
    -config "${CONFIG}" \
    -dataDir "${SNPEFF_DIR}/data" \
    -gff3 \
    -v \
    "${GENOME}"

echo "[snpeff_prep] done — database written to ${DATA_DIR}"
