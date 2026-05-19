#!/usr/bin/env bash
# Run the Snakemake workflow inside the thesis-tools Docker image.
#
# Usage:
#   ./docker/snakemake.sh                        # run default target (bedmethyl)
#   ./docker/snakemake.sh --cores 16             # override core count
#   ./docker/snakemake.sh -n                     # dry-run
#   ./docker/snakemake.sh reads_all --cores 12   # run a specific target
#   ./docker/snakemake.sh variants_all --cores 8 # run clair3 variant calling
#
# All rules including clair3_cpu run inside Docker. Clair3 models are embedded
# in the image at /opt/models/ (copied from hkubal/clair3 during build).
#
# Override image, SnpEff path, or core count via env vars:
#   THESIS_IMAGE=thesis-tools:latest
#   SNPEFF_DIR=/Users/daffa/local/lib/snpEff
#   CORES=8
set -euo pipefail

THESIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNPEFF_DIR="${SNPEFF_DIR:-/Users/daffa/local/lib/snpEff}"
IMAGE="${THESIS_IMAGE:-localhost/thesis-tools:latest}"
CORES="${CORES:-8}"

exec docker run --rm \
    --platform linux/amd64 \
    -v "${THESIS_DIR}:${THESIS_DIR}" \
    -v "${SNPEFF_DIR}:${SNPEFF_DIR}" \
    -w "${THESIS_DIR}/workflow" \
    "${IMAGE}" \
    snakemake --cores "${CORES}" "$@"
