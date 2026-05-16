# Methylation Calling Pipeline Plan
## Dorado Unaligned BAMs → bedMethyl

**Project:** Sorghum bicolor TAA genomic variant study  
**Reference:** GCF_000003195.3 (BTx623, NCBIv3)  
**Sequencing:** ONT R10.4.1 duplex, 400bps SUP, Dorado v1.3.0  

---

## Overview

```
Dorado unaligned BAM (MM/ML tags present)
        │
        ▼
[1] Validate MM/ML tag presence
        │
        ▼
[2] Align to reference (minimap2 + -y flag)
        │
        ▼
[3] Sort & index aligned BAM (samtools)
        │
        ▼
[4] Merge per-library BAMs → per-sample BAM (samtools merge)
        │
        ▼
[5] QC: verify MM/ML tags survived alignment
        │
        ▼
[6] modkit pileup → bedMethyl
        │
        ▼
[7] Filter & summarize bedMethyl output
```

---

## Step 1 — Validate MM/ML Tags in Dorado BAMs

Before aligning, confirm that the Dorado unaligned BAMs carry methylation tags.
This prevents silent data loss downstream.

```bash
# Check MM/ML presence in a sample of reads
samtools view resources/trim_bam/{library}.bam | head -1000 | \
    awk '{found=0; for(i=12;i<=NF;i++) if($i~/^MM:/) found=1; print found}' | \
    sort | uniq -c

# Expected output:
# 1000 1   ← all reads have MM tag (good)
#   0  0   ← no reads lack MM tag (good)

# Also verify the modification type encoded (e.g., 5mC, 5hmC)
samtools view resources/trim_bam/{library}.bam | head -5 | \
    tr '\t' '\n' | grep "^MM:"
# e.g. MM:Z:C+m?,3,1,0;   ← 5mC on cytosines
#      MM:Z:C+h?,2,0,4;   ← 5hmC on cytosines
```

**Expected tags from Dorado SUP duplex with modification calling:**
- `MM:Z:` — modification base positions (base, strand, modification type, skip counts)
- `ML:B:C,` — modification likelihood scores (0–255 per position)

If tags are absent, the Dorado BAMs were basecalled without `--modified-bases` and
must be re-basecalled before proceeding.

---

## Step 2 — Align to Reference (Preserve MM/ML Tags)

The critical step. minimap2's `-y` flag copies all aux tags (MM, ML, mv, etc.)
from the input BAM through to the output aligned BAM.

**Without `-y`, all methylation information is silently discarded.**

```bash
minimap2 \
    -ax map-ont \
    -t {threads} \
    -y \                      # CRITICAL: preserve aux tags (MM, ML, mv)
    --secondary=no \          # suppress secondary alignments (reduces BAM size)
    resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
    resources/trim_bam/{library}.bam \
    | samtools sort -@ {threads} -o resources/align_bam/{library}.bam
```

**Notes:**
- minimap2 accepts unaligned BAM directly as input (no FASTQ conversion needed)
- `--secondary=no` is recommended for methylation calling to avoid double-counting
  modified bases at multi-mapping loci
- Threads: 8–16 recommended per library on the CPU HPC (48 cores available)

---

## Step 3 — Index Aligned BAMs

Required for modkit pileup to perform random-access region queries.

```bash
samtools index resources/align_bam/{library}.bam
```

---

## Step 4 — Merge Per-Library BAMs into Per-Sample BAMs

SBC11 has three libraries (r0075, r0078, r0078-2); others are single-library.
`samtools merge` preserves all aux tags including MM/ML by default.

```bash
# Multi-library sample (SBC11)
samtools merge -f -r -@ {threads} - \
    resources/align_bam/r0075.bam \
    resources/align_bam/r0078.bam \
    resources/align_bam/r0078-2.bam \
    | samtools sort -@ {threads} -o resources/align_bam_sample/SBC11.bam
samtools index resources/align_bam_sample/SBC11.bam

# Single-library samples (SBC04, SBC10, SBC23) — symlink to avoid duplication
ln -sf resources/align_bam/r0074.bam resources/align_bam_sample/SBC04.bam
ln -sf resources/align_bam/r0074.bam.bai resources/align_bam_sample/SBC04.bam.bai
# repeat for SBC10 (r0066) and SBC23 (r0076)
```

**Flag note:** `-r` in `samtools merge` re-attaches read group headers from the
input BAMs, which is useful for tracing reads back to their library of origin.

---

## Step 5 — QC: Verify MM/ML Tags Survived Alignment

Check a sample of reads in the merged BAM before running modkit (which is slow).

```bash
# Quick tag presence check
samtools view resources/align_bam_sample/{sample}.bam | head -1000 | \
    awk '{found=0; for(i=12;i<=NF;i++) if($i~/^MM:/) found=1; print found}' | \
    sort | uniq -c

# Inspect a single read in detail
samtools view resources/align_bam_sample/{sample}.bam | head -1 | tr '\t' '\n' | grep -E "^(MM|ML)"

# Run modkit summary (fast, no alignment required, reports global mod frequency)
modkit summary resources/align_bam_sample/{sample}.bam
```

`modkit summary` output will show per-modification-type read counts and mean
modification frequency. If it shows 0 modified bases, tags are absent.

---

## Step 6 — modkit pileup → bedMethyl

This is the core methylation calling step. modkit aggregates per-read base
modification probabilities into per-genomic-position modification frequencies.

### 6a. Single-sample pileup

```bash
modkit pileup \
    resources/align_bam_sample/{sample}.bam \
    resources/bedmethyl/{sample}.bed \
    --ref resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
    --threads {threads} \
    --preset traditional \        # CpG context (use --cpg for CpG only)
    --filter-threshold 0.75 \     # min modification probability to call (default: 0.75)
    --ignore h \                  # ignore 5hmC if only interested in 5mC
    --log-filepath logs/modkit/{sample}.log
```

**Key flags:**
| Flag | Purpose |
|---|---|
| `--ref` | Reference FASTA; enables reference-anchored CpG grouping |
| `--preset traditional` | CpG context; groups C+m and C+h into combined CpG output |
| `--filter-threshold` | Probability cutoff (0–1); reads below this are counted as "no call" |
| `--ignore h` | Exclude 5hmC from output (relevant if model calls both 5mC and 5hmC) |
| `--threads` | Parallelism; 16–32 recommended on CPU HPC |
| `--cpg` | Restrict output to CpG dinucleotide context only |
| `--combine-strands` | Merge +/- strand CpG counts into single position |

### 6b. Recommended command for this project

For sorghum with 5mC calling (R10.4.1 SUP model), combine strands at CpG sites:

```bash
modkit pileup \
    resources/align_bam_sample/{sample}.bam \
    resources/bedmethyl/{sample}.bed \
    --ref resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
    --threads 24 \
    --cpg \
    --combine-strands \
    --filter-threshold 0.75 \
    --log-filepath logs/modkit/{sample}.log
```

### 6c. bedMethyl output format

Each row in the output represents one genomic position:

```
chrom  start  end  name  score  strand  start  end  color  N_valid  pct_mod  N_mod  N_canonical  N_nocall  N_filtered  N_diff  N_delete
```

| Column | Description |
|---|---|
| chrom | Chromosome (matches BAM reference names) |
| start / end | 0-based genomic position |
| name | Modification type (e.g., `m` for 5mC) |
| score | Same as N_valid (total informative reads) |
| strand | `+` or `-` (or `.` if combined) |
| N_valid | Reads with a confident call (mod or canonical) |
| pct_mod | Percent methylation (0–100) |
| N_mod | Reads called as modified |
| N_canonical | Reads called as canonical |
| N_nocall | Reads below filter threshold (uncertain) |

---

## Step 7 — Filter & Summarize bedMethyl Output

Raw bedMethyl includes positions with very low coverage that are unreliable.
Apply a minimum coverage filter before downstream analysis.

### 7a. Filter by minimum coverage

```bash
# Require at least 10 valid reads per position
awk '$10 >= 10' resources/bedmethyl/{sample}.bed \
    > resources/bedmethyl/{sample}.filtered.bed
```

### 7b. Basic summary statistics

```bash
# Count positions passing filter
wc -l resources/bedmethyl/{sample}.filtered.bed

# Distribution of methylation levels
awk '$10 >= 10 {print $11}' resources/bedmethyl/{sample}.bed | \
    awk 'BEGIN{n=0; s=0} {n++; s+=$1} END{print "Mean pct_mod:", s/n, "N positions:", n}'

# Highly methylated positions (>80%)
awk '$10 >= 10 && $11 > 80' resources/bedmethyl/{sample}.filtered.bed | wc -l
```

### 7c. Cross-sample comparison (optional)

For comparing methylation between high-TAA (SBC10) and low-TAA (SBC11, SBC04) samples:

```bash
# Intersect bedMethyl positions across samples using bedtools
bedtools intersect \
    -a resources/bedmethyl/SBC10.filtered.bed \
    -b resources/bedmethyl/SBC11.filtered.bed \
    -wa -wb \
    > resources/bedmethyl/SBC10_vs_SBC11.intersect.bed
```

---

## Snakemake Integration

The pipeline maps onto two new rules to be added to `reads_preprocessing.smk`,
plus a new `methylation.smk` workflow:

### Rules to add/modify in `reads_preprocessing.smk`

| Rule | Change |
|---|---|
| `align_reads` | Switch input from `.fq` to Dorado `.bam`; add `-y` to minimap2 |
| `symlink_sample_bam` | No change needed |
| `merge_sample_bam` | No change needed (`samtools merge` preserves MM/ML) |

### New `methylation.smk` workflow

```
rule validate_mod_tags      → check MM/ML present in Dorado BAMs
rule modkit_summary         → fast global QC per sample BAM
rule modkit_pileup          → produce per-sample bedMethyl
rule filter_bedmethyl       → apply coverage filter (N_valid ≥ 10)
```

---

## Directory Structure

```
thesis/
├── resources/
│   ├── trim_bam/          # INPUT: Dorado unaligned BAMs (with MM/ML)
│   │   ├── r0066.bam        # SBC10
│   │   ├── r0074.bam        # SBC04
│   │   ├── r0075.bam        # SBC11 lib 1
│   │   ├── r0078.bam        # SBC11 lib 2
│   │   ├── r0078-2.bam      # SBC11 lib 3
│   │   └── r0076.bam        # SBC23
│   ├── align_bam/           # Per-library aligned BAMs (with MM/ML preserved)
│   ├── align_bam_sample/    # Per-sample merged BAMs
│   └── bedmethyl/           # OUTPUT: bedMethyl files
│       ├── SBC04.bed
│       ├── SBC10.bed
│       ├── SBC11.bed
│       └── SBC23.bed
└── workflow/
    ├── reads_preprocessing.smk   # MODIFY: align_reads rule
    └── methylation.smk           # NEW: modkit pileup workflow
```

---

## Key Checkpoints & Common Pitfalls

| Step | Checkpoint | Common Pitfall |
|---|---|---|
| Step 1 | MM/ML tags present in Dorado BAM | Dorado run without `--modified-bases` flag |
| Step 2 | `-y` in minimap2 call | Omitting `-y` silently strips all mod tags |
| Step 2 | minimap2 input is `.bam` not `.fq` | FASTQ cannot carry MM/ML tags |
| Step 5 | `modkit summary` shows non-zero mod counts | Tags dropped during merge or sort |
| Step 6 | `--ref` provided to modkit pileup | Without ref, CpG grouping is disabled |
| Step 7 | Coverage filter applied before analysis | Low-coverage positions inflate variance |

---

## Software Versions

| Tool | Version / Notes |
|---|---|
| Dorado | v1.3.0 (stereo duplex v1.4, SUP model) |
| minimap2 | any recent version; verify `-y` flag supported |
| samtools | ≥1.15 recommended |
| modkit | any version supporting `pileup` and `--cpg` |
| bedtools | for cross-sample intersections |