#!/usr/bin/env bash
set -euo pipefail

VCF_DIR="results/vcf_processing"
OUT_DIR="analysis/data/vcf/private_variants"
LOG_DIR="analysis/logs"
mkdir -p "${LOG_DIR}"
LOG="${LOG_DIR}/private_variants.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

SAMPLES=(SBC4 SBC10 SBC11 SBC23)
VCFS=()
for s in "${SAMPLES[@]}"; do
    VCFS+=("${VCF_DIR}/${s}.phased.vcf.gz")
done

# Index any VCF that is missing a .tbi
for vcf in "${VCFS[@]}"; do
    if [[ ! -f "${vcf}.tbi" ]]; then
        echo "Indexing ${vcf}..."
        bcftools index -t "$vcf"
    fi
done

mkdir -p "$OUT_DIR"

echo "Running bcftools isec for private variants..."
bcftools isec -n =1 -p "$OUT_DIR" "${VCFS[@]}"

# Rename outputs to sample names for clarity
for i in "${!SAMPLES[@]}"; do
    idx=$(printf "%04d" "$i")
    src="${OUT_DIR}/${idx}.vcf"
    dst="${OUT_DIR}/${SAMPLES[$i]}.private.vcf"
    if [[ -f "$src" ]]; then
        mv "$src" "$dst"
        bgzip -f "$dst"
        bcftools index -t "${dst}.gz"
        echo "${SAMPLES[$i]}: $(bcftools view -H ${dst}.gz | wc -l) private variants"
    fi
done

echo "Done. Output in ${OUT_DIR}/"
