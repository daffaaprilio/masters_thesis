#!/usr/bin/env bash
# Assign CpG / CHG / CHH sequence context to 5mC bedmethyl positions using
# bedtools getfasta, then aggregate counts per sample and context.
#
# Output TSV columns:
#   sample  context  n_sites  mean_pct_mod  total_N_mod  total_N_valid
#
# Usage:
#   ./docker/run.sh bash analysis/scripts/methylation_context_samples.sh

set -euo pipefail

REF="resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
BED_DIR="resources/bedmethyl"
OUT_DIR="analysis/data/methylation_context"
LOG_DIR="analysis/logs"
MIN_COV=10
SAMPLES=(SBC4 SBC10 SBC11 SBC23)

mkdir -p "$OUT_DIR" "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/methylation_context_${TIMESTAMP}.log"
OUT_TSV="$OUT_DIR/methylation_context_${TIMESTAMP}.tsv"

log() { echo "$1" | tee -a "$LOG_FILE"; }

log "Methylation context assignment (bedtools getfasta)"
log "Started:     $(date '+%Y-%m-%d %H:%M:%S')"
log "Reference:   $REF"
log "Min coverage: $MIN_COV"
log ""

printf "sample\tcontext\tn_sites\tmean_pct_mod\ttotal_N_mod\ttotal_N_valid\n" > "$OUT_TSV"

for sample in "${SAMPLES[@]}"; do
    log "════════════════════════════════════════════════════════════"
    log "  $sample"

    # Prefer filtered bedmethyl, fall back to raw
    bed=""
    for fname in "${BED_DIR}/${sample}.filtered.bed" "${BED_DIR}/${sample}.bed"; do
        [[ -f "$fname" ]] && { bed="$fname"; break; }
    done

    if [[ -z "$bed" ]]; then
        log "  WARNING: no bedmethyl file found — skipping"
        continue
    fi
    log "  Input: $bed"

    tmp_filtered=$(mktemp --suffix=.bed)
    tmp_windows=$(mktemp --suffix=.bed)
    tmp_fasta=$(mktemp --suffix=.fa)

    # ── 1. Filter for 5mC rows with sufficient coverage ──────────────────────
    # bedmethyl columns: chrom(1) start(2) end(3) name(4) score(5) strand(6)
    #   thick_start(7) thick_end(8) color(9) N_valid(10) pct_mod(11) N_mod(12) ...
    awk -v mincov="$MIN_COV" '$4 == "m" && $10 >= mincov' "$bed" > "$tmp_filtered"
    n_sites=$(wc -l < "$tmp_filtered")
    log "  5mC sites (N_valid ≥ $MIN_COV): $(printf "%'d" "$n_sites")"

    # ── 2. Build strand-aware 3-base windows for bedtools getfasta ───────────
    # + strand C at pos i: extract [i, i+3)  → seq is C + next_2_bases
    # - strand C at pos i: extract [i-2, i+1) with -s → reverse complement
    #   gives C + complement(ref[i-1]) + complement(ref[i-2])
    awk 'BEGIN{OFS="\t"} {
        chrom = $1; start = $2; strand = $6
        if (strand == "+") {
            wstart = start
            wend   = start + 3
        } else {
            wstart = (start - 2 >= 0 ? start - 2 : 0)
            wend   = start + 1
        }
        print chrom, wstart, wend, NR, ".", strand
    }' "$tmp_filtered" > "$tmp_windows"

    # ── 3. Extract sequences (strand-aware reverse complement on - strand) ───
    bedtools getfasta -fi "$REF" -bed "$tmp_windows" -s -fo "$tmp_fasta" 2>>"$LOG_FILE"

    # ── 4. Assign CpG / CHG / CHH and aggregate ──────────────────────────────
    # Parse sequences: b1=seq[2], b2=seq[3]
    #   b1 == G → CpG
    #   b2 == G → CHG
    #   else    → CHH
    # Sequences shorter than 3 bases (chromosome edge) are classified CHH.
    paste \
        <(awk 'BEGIN{OFS="\t"} {print $10, $11, $12}' "$tmp_filtered") \
        <(grep -v "^>" "$tmp_fasta" | awk '{
            seq = toupper($0)
            b1  = substr(seq, 2, 1)
            b2  = substr(seq, 3, 1)
            if      (b1 == "G") print "CpG"
            else if (b2 == "G") print "CHG"
            else                print "CHH"
        }') \
    | awk -v sample="$sample" '
        BEGIN { OFS = "\t" }
        {
            ncov = $1; pct = $2; nmod = $3; ctx = $4
            n[ctx]++
            sum_pct[ctx]  += pct
            sum_nmod[ctx] += nmod
            sum_ncov[ctx] += ncov
        }
        END {
            for (ctx in n)
                printf "%s\t%s\t%d\t%.4f\t%d\t%d\n",
                    sample, ctx, n[ctx],
                    sum_pct[ctx] / n[ctx],
                    sum_nmod[ctx], sum_ncov[ctx]
        }
    ' >> "$OUT_TSV"

    rm -f "$tmp_filtered" "$tmp_windows" "$tmp_fasta"
    log "  Done"
done

log ""
log "Output TSV: $OUT_TSV"
log "Log:        $LOG_FILE"
log "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
