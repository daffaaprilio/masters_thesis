#!/usr/bin/env bash
# Run any bioinformatics command through the thesis-tools Docker image.
# Mounts the thesis directory and the local SnpEff data directory.
#
# Usage:
#   ./docker/run.sh samtools view -H sample.bam
#   ./docker/run.sh modkit pileup sample.bam out.bed --ref ref.fna --threads 8
#   ./docker/run.sh snpEff ann -config /path/snpEff.config Sorghum_bicolor in.vcf
#
# Override image or SnpEff path via env vars:
#   THESIS_IMAGE=thesis-tools:latest
#   SNPEFF_DIR=/custom/path/to/snpEff
set -euo pipefail

THESIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNPEFF_DIR="${SNPEFF_DIR:-${THESIS_DIR}/resources/snpeff}"
IMAGE="${THESIS_IMAGE:-localhost/thesis-tools:latest}"

exec docker run --rm \
    --platform linux/amd64 \
    --security-opt label=disable \
    -v "${THESIS_DIR}:${THESIS_DIR}" \
    -w "$(pwd)" \
    "${IMAGE}" \
    "$@"
