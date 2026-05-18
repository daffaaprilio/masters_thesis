#!/usr/bin/env bash
# Build the thesis-tools Docker image.
# Run from anywhere; always uses the thesis root as build context.
set -euo pipefail

THESIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

docker build \
    --platform linux/amd64 \
    -f "${THESIS_DIR}/docker/Dockerfile" \
    -t thesis-tools:latest \
    "${THESIS_DIR}"

echo "Built thesis-tools:latest"
