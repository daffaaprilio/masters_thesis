#!/bin/zsh

WDIR="$(cd "$(dirname "$0")/.." && pwd)"

SNPEFF_DIR="/Users/daffa/local/lib/snpEff"
SNPEFF_DB="Sorghum_bicolor"
RENAME_MAP="${WDIR}/workflow/scripts/synonyms.txt"

OUT_DIR="${WDIR}/analysis/data/vcf/annotated"

mkdir -p "${OUT_DIR}" "${WDIR}/analysis/logs"

SAMPLE=$1
INPUT_VCF=$2
LOG="${WDIR}/analysis/logs/annotate.${SAMPLE}.$(date +%Y%m%d_%H%M%S).log"

(
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
) > "${LOG}" 2>&1
