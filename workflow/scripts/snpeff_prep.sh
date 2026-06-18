#!/usr/bin/env bash
# Build a custom SnpEff database for Sorghum bicolor from the same NCBI GTF
# used by SIFT4G, ensuring consistent gene models between the two tools.
#
# Usage (from any directory within the repo):
#   bash workflow/scripts/snpeff_prep.sh
set -euo pipefail

THESIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SNPEFF_DIR="${THESIS_DIR}/resources/snpeff"

DB="Sorghum_bicolor_NCBIv3"
GTF="${THESIS_DIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gtf"
REF="${THESIS_DIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"

DATA_DIR="${SNPEFF_DIR}/data"
DB_DIR="${DATA_DIR}/${DB}"
CONFIG="${SNPEFF_DIR}/snpEff.config"

mkdir -p "${DB_DIR}"

echo "[snpeff_prep] copying GTF annotation..."
cp "${GTF}" "${DB_DIR}/genes.gtf"

echo "[snpeff_prep] linking genome FASTA..."
ln -sf "${REF}" "${DB_DIR}/sequences.fa"

echo "[snpeff_prep] writing snpEff.config..."
cat > "${CONFIG}" <<EOF
data.dir = ${DATA_DIR}

${DB}.genome : Sorghum bicolor NCBIv3
EOF

echo "[snpeff_prep] building SnpEff database from GTF..."
snpEff build \
    -gtf22 \
    -v \
    -noCheckCds \
    -noCheckProtein \
    -config "${CONFIG}" \
    "${DB}"

echo "[snpeff_prep] done — database written to ${DB_DIR}"
