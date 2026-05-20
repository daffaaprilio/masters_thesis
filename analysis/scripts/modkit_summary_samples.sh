#!/usr/bin/env bash
# Run modkit summary on per-sample BAMs and collect methylation status.
# Outputs a TSV summary to analysis/data/modkit_summary/ and a timestamped log.
#
# Usage:
#   ./docker/run.sh bash analysis/scripts/modkit_summary_samples.sh
#   # or directly (if modkit is on PATH):
#   bash analysis/scripts/modkit_summary_samples.sh

set -euo pipefail

BAM_DIR="resources/align_bam_sample"
OUT_DIR="analysis/data/modkit_summary"
LOG_DIR="analysis/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/modkit_summary_${TIMESTAMP}.log"
TSV_FILE="$OUT_DIR/modkit_summary_${TIMESTAMP}.tsv"

log() { echo "$1" | tee -a "$LOG_FILE"; }

log "modkit summary — per-sample BAMs"
log "Started: $(date '+%Y-%m-%d %H:%M:%S')"
log "BAM dir: $BAM_DIR"
log ""

# Write TSV header
printf "sample\tbase\tcode\tpass_count\tpass_frac\tall_count\tall_frac\n" \
    > "$TSV_FILE"

for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    log "════════════════════════════════════════════════════════════"
    log "  $sample"
    log "════════════════════════════════════════════════════════════"

    # Capture modkit summary output; it prints a table to stdout
    summary_output=$(modkit summary "$bam" 2>>"$LOG_FILE")
    log "$summary_output"
    log ""

    # Parse data rows: first column is base (A/C uppercase), second is mod code.
    # modkit summary columns: base  code  pass_count  pass_frac  all_count  all_frac
    echo "$summary_output" | awk -v sample="$sample" '
        /^[[:space:]]+[AC][[:space:]]/ {
            gsub(/^[[:space:]]+/, "")
            n = split($0, f, /[[:space:]]+/)
            if (n >= 6) {
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    sample, f[1], f[2], f[3], f[4], f[5], f[6]
            }
        }
    ' >> "$TSV_FILE"
done

log "════════════════════════════════════════════════════════════"
log "Summary TSV: $TSV_FILE"
log "Log:         $LOG_FILE"
log "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
