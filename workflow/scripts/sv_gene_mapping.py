#!/opt/conda/bin/python3
"""
Map structural variants to the gene(s) they affect — the SV counterpart of
genomics_scoring.py's front half, but WITHOUT scoring.

Gene assignment comes from the SnpEff ANN field (the SV group VCF must already be
SnpEff-annotated; see the annotate_sv rule), not from raw GFF interval overlap.
SnpEff classifies each SV's consequence per gene — transcript_ablation, exon_loss,
frameshift, feature_fusion (BND), duplication / inversion (DUP/INV), plus
upstream/downstream/intron MODIFIERs — and emits one ANN entry per affected
transcript. We reduce those to ONE entry per gene (most severe impact wins) and
emit long format: ONE ROW PER SV–GENE. A deletion hitting five genes → five rows
(same sv_id).

Genes are resolved against the NCBI gene_info Symbol set (same approach as
build_ranked_genes.py): ANN gene names that are real gene symbols are kept;
chromosome-level summary entries (empty gene) and compound intergenic names
(e.g. "LOCa-LOCb") are dropped, so an SV with no affected single gene yields no
rows. The effect/impact are recorded as descriptive columns — no numeric
score/rank, the deliberate difference from genomics_scoring.py.

Input is one SnpEff-annotated SV group VCF (results/snpeff_sv/{group}.annotated.vcf.gz);
the group label is derived from the filename, mirroring genomics_scoring.py.
"""

import argparse
import os
import sys

import cyvcf2
import pandas as pd

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]  # positional, == combined VCF columns

# SnpEff ANN sub-fields (pipe-delimited), same layout as genomics_scoring.ANN_FIELDS.
ANN_GENE_NAME    = 3   # Gene_Name  (LOC* for real genes)
ANN_EFFECT       = 1   # Annotation (effect)
ANN_IMPACT       = 2   # Annotation_Impact
ANN_FEATURE_TYPE = 5   # Feature_Type
ANN_BIOTYPE      = 7   # Transcript_BioType

# Impact severity for reducing many ANN entries per gene to one (HIGH = most severe).
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

OUTPUT_COLS = [
    "group", "chrom", "pos", "end", "sv_id", "svtype", "svlen", "support", "af",
    *(f"gt_{s}" for s in SAMPLES),
    "gene_id", "effect", "impact", "feature_type", "biotype",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="SnpEff-annotated SV group VCF (results/snpeff_sv/<group>.annotated.vcf.gz)",
    )
    parser.add_argument(
        "--gene-info", required=True,
        help="NCBI gene_info tab file (resources/NCBI_FTP/gene_info_4558)",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory (TSV named <group>.sv_genes.tsv)",
    )
    return parser.parse_args()


def load_valid_symbols(gene_info_path: str) -> set:
    """Set of real gene symbols (LOC*, tRNA names, …) from NCBI gene_info.

    Mirrors build_ranked_genes.load_sorbi_to_loc: a SnpEff ANN gene name is kept
    only if it is one of these, which rejects compound intergenic names
    ("LOCa-LOCb") and empty chromosome-level entries.
    """
    gi = pd.read_csv(gene_info_path, sep="\t", usecols=["Symbol"])
    return set(gi["Symbol"].dropna().astype(str))


def reduce_to_one_per_gene(ann_raw: str, valid_symbols: set) -> dict:
    """Collapse a record's ANN entries to ONE per gene: most severe impact wins.

    Returns {gene_id: (effect, impact, feature_type, biotype)} for genes that are
    real symbols. Ties broken by ANN order (first/most-severe SnpEff entry).
    """
    best: dict[str, tuple] = {}
    for order, entry in enumerate(ann_raw.split(",")):
        parts = entry.split("|")
        if len(parts) <= ANN_BIOTYPE:
            continue
        gene = parts[ANN_GENE_NAME]
        if gene not in valid_symbols:
            continue
        impact = parts[ANN_IMPACT]
        rank = (IMPACT_RANK.get(impact, 99), order)
        if gene not in best or rank < best[gene][0]:
            best[gene] = (
                rank,
                parts[ANN_EFFECT],
                impact,
                parts[ANN_FEATURE_TYPE],
                parts[ANN_BIOTYPE],
            )
    return {g: v[1:] for g, v in best.items()}


def map_svs_to_genes(vcf_path: str, valid_symbols: set, group: str):
    """Yield one dict per SV–gene pair; returns (rows, n_sv, n_with_gene)."""
    rows = []
    n_sv = 0
    n_with_gene = 0

    vcf = cyvcf2.VCF(vcf_path)
    sample_idx = {s: vcf.samples.index(s) for s in SAMPLES if s in vcf.samples}

    for v in vcf:
        n_sv += 1

        ann_raw = v.INFO.get("ANN")
        if not ann_raw:
            continue
        gene_hits = reduce_to_one_per_gene(str(ann_raw), valid_symbols)
        if not gene_hits:
            continue
        n_with_gene += 1

        end_raw = v.INFO.get("END")
        end = int(end_raw) if end_raw is not None else v.POS

        gts = {}
        for s in SAMPLES:
            j = sample_idx.get(s)
            if j is None:
                gts[s] = "./."
                continue
            a = v.genotypes[j]
            sep = "|" if a[2] else "/"
            gts[s] = f"{a[0]}{sep}{a[1]}"

        base = {
            "group":   group,
            "chrom":   v.CHROM,
            "pos":     v.POS,
            "end":     end,
            "sv_id":   v.ID,
            "svtype":  v.INFO.get("SVTYPE"),
            "svlen":   v.INFO.get("SVLEN"),
            "support": v.INFO.get("SUPPORT"),
            "af":      v.INFO.get("AF"),
            **{f"gt_{s}": gts[s] for s in SAMPLES},
        }

        for gene_id, (effect, impact, feature_type, biotype) in gene_hits.items():
            rows.append({
                **base,
                "gene_id":      gene_id,
                "effect":       effect,
                "impact":       impact,
                "feature_type": feature_type,
                "biotype":      biotype,
            })

    return rows, n_sv, n_with_gene


def main():
    args = parse_args()

    group = os.path.basename(args.input).split(".")[0]
    os.makedirs(args.output, exist_ok=True)

    print(f"[sv_gene_mapping] Loading gene symbols: {args.gene_info}", file=sys.stderr)
    valid_symbols = load_valid_symbols(args.gene_info)
    print(f"[sv_gene_mapping]   {len(valid_symbols):,} valid gene symbols", file=sys.stderr)

    print(f"[sv_gene_mapping] Mapping SVs (SnpEff ANN): {args.input} (group={group})",
          file=sys.stderr)
    rows, n_sv, n_with_gene = map_svs_to_genes(args.input, valid_symbols, group)

    df = pd.DataFrame(rows, columns=OUTPUT_COLS)
    out_tsv = os.path.join(args.output, f"{group}.sv_genes.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"[sv_gene_mapping]   SVs read              : {n_sv:,}", file=sys.stderr)
    print(f"[sv_gene_mapping]   SVs with a gene hit   : {n_with_gene:,}", file=sys.stderr)
    print(f"[sv_gene_mapping]   SVs without gene hit  : {n_sv - n_with_gene:,}", file=sys.stderr)
    print(f"[sv_gene_mapping]   SV–gene rows          : {len(df):,}", file=sys.stderr)
    print(f"[sv_gene_mapping] Written → {out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
