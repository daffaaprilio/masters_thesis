#!/bin/zsh

WDIR="$(cd "$(dirname "$0")/.." && pwd)"

IN_DIR="${WDIR}/results/vcf_processing"
OUT_DIR="${WDIR}/analysis/data/vcf/merged"
LOG="${WDIR}/analysis/logs/merge_samples.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "${OUT_DIR}" "${WDIR}/analysis/logs"

SAMPLES=(SBC10 SBC11 SBC4 SBC23)
INPUT_VCFS=("${(@)SAMPLES/#/${IN_DIR}/}" )
INPUT_VCFS=("${(@)INPUT_VCFS/%/.phased.vcf.gz}")

(
    echo "[$(date +%T)] Merging ${#SAMPLES[@]} samples: ${SAMPLES[*]}..."
    bcftools merge \
        -O z \
        -o "${OUT_DIR}/merged.phased.vcf.gz" \
        "${INPUT_VCFS[@]}"

    bcftools index "${OUT_DIR}/merged.phased.vcf.gz"
    echo "[$(date +%T)] Done: ${OUT_DIR}/merged.phased.vcf.gz"
) > "${LOG}" 2>&1
