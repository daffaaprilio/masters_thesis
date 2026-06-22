#!/usr/bin/env bash
# TEMPORARY orchestrator — variant-group annotation chain WITHOUT Snakemake.
# Run inside the toolchain container:  ./docker/run.sh bash run_variant_groups.sh
# For all 15 sample-sharing groups: bcftools isec -> SnpEff -> SIFT4G, then UpSet plot.
# Idempotent: skips any per-group step whose output already exists.
set -euo pipefail

WDIR="$(cd "$(dirname "$0")" && pwd)"; cd "$WDIR"
SAMPLES=(SBC4 SBC10 SBC11 SBC23)
PHASED_DIR=results/vcf_processing
VARGROUP_DIR=results/variant_groups
SNPEFF_DIR=results/snpeff
SIFT_DIR=results/sift4g
SIFT_DB=resources/sift4g/sorghum_db/NCBIv3
SIFT_JAR=/opt/SIFT4G_Annotator.jar
log() { echo "[$(date +%T)] $*"; }

# ── Prerequisites (fail fast) ──────────────────────────────────────────────────
for s in "${SAMPLES[@]}"; do
  [ -f "$PHASED_DIR/$s.phased.vcf.gz" ] || { echo "missing phased VCF: $s"; exit 1; }
done
[ -f resources/snpeff/snpEff.config ]        || { echo "SnpEff DB not built";  exit 1; }
[ -f resources/sift4g/sorghum_db/.db_built ] || { echo "SIFT4G DB not built";  exit 1; }
[ -f "$SIFT_JAR" ]                           || { echo "SIFT jar missing";     exit 1; }
mkdir -p "$VARGROUP_DIR" "$SNPEFF_DIR" "$SIFT_DIR"
VCFS=(); for s in "${SAMPLES[@]}"; do VCFS+=("$PHASED_DIR/$s.phased.vcf.gz"); done

# ── Per-group chain: isec -> SnpEff -> SIFT4G ──────────────────────────────────
for mask in $(seq 1 15); do
  pattern=""; label=""; w=0
  for i in 0 1 2 3; do
    if [ $(( (mask >> i) & 1 )) -eq 1 ]; then   # [ $((..)) -eq 1 ], not bare (( )) — set -e safe
      pattern+="1"
      if [ -z "$label" ]; then label="${SAMPLES[$i]}"; w=$((i+1)); else label="${label}_${SAMPLES[$i]}"; fi
    else
      pattern+="0"
    fi
  done
  log "=== $label (pattern=$pattern -w$w) ==="
  gvcf="$VARGROUP_DIR/$label.vcf.gz"
  avcf="$SNPEFF_DIR/$label.annotated.vcf.gz"
  svcf="$SIFT_DIR/$label.sift4g.vcf.gz"

  # 1) intersection
  if [ -f "$gvcf" ]; then log "  isec   skip"; else
    log "  isec"
    bcftools isec -n ~"$pattern" -w"$w" -O z -o "$gvcf" "${VCFS[@]}"
    bcftools index -t "$gvcf"
  fi

  # 2) SnpEff (reuses annotate_vcf.sh -> {label}.annotated.vcf.gz + .csi)
  if [ -f "$avcf" ]; then log "  snpeff skip"; else
    log "  snpeff"
    bash workflow/scripts/annotate_vcf.sh "$label" "$gvcf" "$SNPEFF_DIR"
    bcftools view -h "$avcf" >/dev/null   # sanity: snpEff|bgzip in that script has no pipefail
  fi

  # 3) SIFT4G (inlined from rule annotate_sift)
  if [ -f "$svcf" ]; then log "  sift4g skip"; else
    log "  sift4g"
    tmp="$SIFT_DIR/${label}_tmp"; stem="$label.annotated"; mkdir -p "$tmp"
    bcftools view "$avcf" -O v -o "$tmp/$stem.vcf"
    java -jar "$SIFT_JAR" -c -i "$tmp/$stem.vcf" -d "$SIFT_DB" -r "$SIFT_DIR"
    sift_vcf=$(ls "$SIFT_DIR/${stem}_SIFTpredictions"*.vcf 2>/dev/null | head -1 || true)
    [ -z "${sift_vcf:-}" ] && sift_vcf=$(ls "$SIFT_DIR/${stem}"*SIFT*.vcf 2>/dev/null | head -1 || true)
    [ -n "${sift_vcf:-}" ] || { echo "SIFT produced no VCF for $label"; exit 1; }
    bcftools sort "$sift_vcf" | bgzip -c > "$svcf"
    bcftools index -t "$svcf"
    rm -f "$sift_vcf"; rm -rf "$tmp"
  fi
done

# ── UpSet plot from the group VCF counts (always regenerate; cheap) ────────────
log "=== UpSet plot ==="
python3 workflow/scripts/variant_upset.py

log "DONE — 15 groups annotated; figure at analysis/figures/variant_upset.png"
