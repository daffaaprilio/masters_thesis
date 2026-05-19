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
printf "sample\tbase\tstrand\tn_called\tn_mod\tpct_modified\tn_canonical\tn_other\tn_delete\tn_fail\tn_diff\tn_nocall\n" \
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

    # Parse the table lines (skip header and separator lines)
    # modkit summary table columns (as of modkit 0.2.x):
    #   base  strand  n_called  n_mod  pct_modified  n_canonical  n_other_mod  n_delete  n_fail  n_diff  n_no_call
    echo "$summary_output" | awk -v sample="$sample" '
        /^[[:space:]]*[mah]/ {
            # Strip leading whitespace, split on whitespace
            $0 = $0
            gsub(/^[[:space:]]+/, "")
            n = split($0, f, /[[:space:]]+/)
            if (n >= 7) {
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    sample, f[1], f[2], f[3], f[4], f[5], f[6], f[7],
                    (n>=8 ? f[8] : "NA"), (n>=9 ? f[9] : "NA"),
                    (n>=10 ? f[10] : "NA"), (n>=11 ? f[11] : "NA")
            }
        }
    ' >> "$TSV_FILE"
done

log "════════════════════════════════════════════════════════════"
log "Summary TSV: $TSV_FILE"
log "Log:         $LOG_FILE"
log "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
