# Filtering Strategy for Variants Private to SBC10

## Context

I have a SnpEff-annotated VCF file of variants private to SBC10, a high tannin accession of
*Sorghum bicolor*. The goal is to filter this file to retain only high-confidence candidate
variants that may explain elevated trans-aconitic acid (TAA) production in SBC10.

**Input file:** `/home/daffa/Work/2026/thesis/SBC10_private.ann.vcf.gz`  
**Output file:** `/home/daffa/Work/2026/thesis/SBC10_candidate_variants.tsv`

---

## Biological Rationale

The variants of interest are those occurring within candidate genes encoding enzymes that
operate in the TAA biosynthetic pathway. This defines two key properties: the variant must
fall within a protein-coding transcript, and it must be private to SBC10. Additional
properties further narrow the list to biologically credible, technically reliable calls.

---

## Filter Criteria

Retain a variant **only if ALL of the following are satisfied**:

### 1. Variant Call Quality
- `QUAL ≥ 20`

### 2. Read Depth in SBC10
- `10 ≤ DP ≤ 100`
- DP is taken from the SBC10 sample column FORMAT field

### 3. Transcript Biotype
- Must be `protein_coding`
- Sourced from SnpEff ANN field, **subfield index 7** (0-based)
- Rationale: TAA biosynthesis is an enzymatic pathway; only protein-coding transcripts
  produce translated gene products with interpretable functional consequences

### 4. Functional Impact
- Must be `HIGH` or `MODERATE`
- Sourced from SnpEff ANN field, **subfield index 2** (0-based)
- Note: this is a **heuristic classification** by SnpEff reflecting predicted transcript
  disruption, not experimentally validated functional effect
  - `HIGH` — severely disruptive: frameshift, stop gained/lost, splice donor/acceptor loss
  - `MODERATE` — potentially meaningful: missense variant, in-frame indel

### 5. Genotype in SBC10
- Must be homozygous ALT: `GT = 1/1`
- Rationale: SBC10 shows consistently elevated TAA; a fully fixed coding change is a more
  parsimonious explanation of a strong phenotype than a heterozygous call where one
  functional copy may still be present

### 6. GT/AF Concordance in SBC10
- For a `1/1` call, the within-sample ALT allele frequency must be `AF ≥ 0.75`
- AF is taken from the SBC10 sample column FORMAT field
- Rationale: AF is computed directly from raw read counts; a homozygous GT call with a
  low AF indicates a potentially miscalled or technically unreliable variant

---

## ANN Field Structure

The SnpEff ANN field is pipe-delimited (`|`) with 16 subfields per annotation entry.
Multiple transcript annotations per variant are comma-separated.

| Index (0-based) | Field |
|---|---|
| 0 | Allele |
| 1 | Effect (e.g. `missense_variant`) |
| 2 | Impact (`HIGH / MODERATE / LOW / MODIFIER`) |
| 3 | Gene name |
| 4 | Gene ID |
| 5 | Feature type |
| 6 | Feature ID (transcript) |
| 7 | Transcript biotype |
| 8 | Rank |
| 9 | HGVS.c |
| 10 | HGVS.p |
| 11 | cDNA position |
| 12 | CDS position |
| 13 | Protein position |
| 14 | Distance to feature |
| 15 | Errors/warnings |

### Handling Multiple Transcript Annotations

When multiple transcripts are annotated per variant (comma-separated ANN entries):

1. Select the annotation with the **highest impact tier** (HIGH > MODERATE > LOW > MODIFIER)
2. If impact is equal across entries, retain the **first entry**

---

## Output Format

A TSV file with the following columns, in this order:

| Column | Source |
|---|---|
| `chrom` | CHROM field |
| `pos` | POS field |
| `ref` | REF field |
| `alt` | ALT field |
| `qual` | QUAL field |
| `gene` | ANN subfield index 3 |
| `effect` | ANN subfield index 1 |
| `impact` | ANN subfield index 2 |
| `biotype` | ANN subfield index 7 |
| `hgvs_c` | ANN subfield index 9 |
| `hgvs_p` | ANN subfield index 10 |
| `GT` | SBC10 FORMAT: GT |
| `GQ` | SBC10 FORMAT: GQ |
| `DP` | SBC10 FORMAT: DP |
| `AF` | SBC10 FORMAT: AF |

### Sort order
1. Impact tier — HIGH first, then MODERATE
2. Gene name — alphabetical
3. Position — ascending

---

## Notes

- The input VCF is already SBC10-private (produced by `bcftools isec`), so privateness is
  guaranteed upstream and does not need to be re-checked from other sample columns.
- The ANN field is located in the INFO column, prefixed with `ANN=`.
- Print a summary count of passing variants to stdout before saving the TSV.
