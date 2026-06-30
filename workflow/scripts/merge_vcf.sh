#!/usr/bin/env bash

WDIR="$(cd "$(dirname "$0")/../.." && pwd)"

IN_DIR="${WDIR}/results/vcf_processing"
OUT_DIR="${WDIR}/results/vcf_processing"
REF="${WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"

mkdir -p "${OUT_DIR}"

# Canonical SAMPLES order so the merged VCF's sample columns (GT[0..3]) stay stable
# for any downstream genotype-based filtering.
SAMPLES=(SBC4 SBC10 SBC11 SBC23)
INPUT_VCFS=()
for s in "${SAMPLES[@]}"; do
    INPUT_VCFS+=("${IN_DIR}/${s}.phased.vcf.gz")
done

echo "[$(date +%T)] Merging ${#SAMPLES[@]} samples: ${SAMPLES[*]}..."
# Split multiallelics (and left-align indels) right after the merge. This is REQUIRED,
# not cosmetic: bcftools merge collapses different per-sample ALT alleles at one POS
# into a single multiallelic record, so a downstream GT="alt" test would read that site
# as "present" in every sample carrying *any* alt — conflating distinct alleles and
# mis-assigning the sample-sharing groups. Splitting restores per-allele membership that
# matches `bcftools isec` (validated to within indel-realignment noise, ~0.015%).
bcftools merge -O u "${INPUT_VCFS[@]}" \
    | bcftools norm -m -any -f "${REF}" -O z -o "${OUT_DIR}/merged.phased.vcf.gz"

bcftools index -t "${OUT_DIR}/merged.phased.vcf.gz"
echo "[$(date +%T)] Done: ${OUT_DIR}/merged.phased.vcf.gz"
