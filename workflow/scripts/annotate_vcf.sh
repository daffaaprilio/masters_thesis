#!/usr/bin/env bash

WDIR="$(cd "$(dirname "$0")/../.." && pwd)"

SNPEFF_DIR="${WDIR}/resources/snpeff"
SNPEFF_DB="Sorghum_bicolor"
RENAME_MAP="${WDIR}/workflow/scripts/synonyms.txt"

SAMPLE=$1
INPUT_VCF=$2
OUT_DIR="${3:-${WDIR}/analysis/data/vcf/annotated}"

mkdir -p "${OUT_DIR}"

echo "[$(date +%T)] Renaming chromosomes and annotating ${SAMPLE}..."
bcftools annotate \
    --rename-chrs "${RENAME_MAP}" \
    -O v "${INPUT_VCF}" \
    | snpEff ann \
        -config  "${SNPEFF_DIR}/snpEff.config" \
        -dataDir "${SNPEFF_DIR}/data" \
        -v \
        -nodownload \
        -stats   "${OUT_DIR}/${SAMPLE}.stats.html" \
        -csvStats "${OUT_DIR}/${SAMPLE}.stats.csv" \
        "${SNPEFF_DB}" \
    | bgzip -c > "${OUT_DIR}/${SAMPLE}.annotated.vcf.gz"

bcftools index "${OUT_DIR}/${SAMPLE}.annotated.vcf.gz"
echo "[$(date +%T)] Done: ${SAMPLE}"
