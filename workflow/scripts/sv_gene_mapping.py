#!/opt/conda/bin/python3
"""
Map structural variants to the gene(s) they overlap — the SV counterpart of
genomics_scoring.py's front half, but WITHOUT scoring.

The SV VCFs in results/sv_groups/ are vanilla (no SnpEff ANN field), so genes are
assigned by direct interval overlap against the NCBIv3 GFF rather than by parsing
an ANN field. Each SV's span [POS, END] is overlapped with every gene's genomic
coordinates; INS/BND collapse to a point at POS. The GFF seqids (NC_*/NW_*) match
the SV VCF contigs exactly and genes carry Name=LOC*, so the emitted gene IDs line
up 1:1 with the LOC* IDs used elsewhere in the pipeline.

Output is long format: ONE ROW PER SV–GENE OVERLAP. A deletion spanning five
genes produces five rows (same sv_id). SVs overlapping no gene are dropped
(counted to stderr). No score/rank/percentile columns — that is the deliberate
difference from genomics_scoring.py.

Input is one per-group VCF (e.g. results/sv_groups/SBC10.vcf.gz); the group label
is derived from the filename, mirroring genomics_scoring.py.
"""

import argparse
import gzip
import os
import re
import sys

import cyvcf2
import numpy as np
import pandas as pd

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]  # positional, == combined VCF columns

OUTPUT_COLS = [
    "group", "chrom", "pos", "end", "sv_id", "svtype", "svlen", "support", "af",
    *(f"gt_{s}" for s in SAMPLES),
    "gene_id", "gene_biotype", "gene_start", "gene_end", "overlap_bp",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Per-group SV VCF (e.g. results/sv_groups/SBC10.vcf.gz)",
    )
    parser.add_argument(
        "-g", "--gff", required=True,
        help="NCBIv3 GFF3 for gene coordinates",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory (TSV named <group>.sv_genes.tsv)",
    )
    return parser.parse_args()


def load_genes(gff_path: str) -> dict:
    """Per-chromosome gene index for vectorized interval overlap.

    Parses 'gene' and 'pseudogene' features (mirrors genomics_scoring.load_gene_length),
    keeping the genomic span instead of just the length. Keyed by chromosome:

        {chrom: {"start": np.int64[],  # 1-based inclusive (GFF col 4)
                 "end":   np.int64[],  # 1-based inclusive (GFF col 5)
                 "name":  np.str_[],   # Name= attribute (LOC*)
                 "biotype": np.str_[]}}  # gene_biotype= attribute

    Within each chromosome the arrays are sorted by start so downstream output is
    naturally ordered. tRNA/combo keys never appear here, so SVs overlapping only
    those features fall through as intergenic.
    """
    starts:   dict[str, list] = {}
    ends:     dict[str, list] = {}
    names:    dict[str, list] = {}
    biotypes: dict[str, list] = {}

    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = re.search(r"Name=([^;]+)", f[8])
            if not m:
                continue
            chrom = f[0]
            b = re.search(r"gene_biotype=([^;]+)", f[8])
            starts.setdefault(chrom, []).append(int(f[3]))
            ends.setdefault(chrom, []).append(int(f[4]))
            names.setdefault(chrom, []).append(m.group(1))
            biotypes.setdefault(chrom, []).append(b.group(1) if b else "")

    index = {}
    for chrom in starts:
        s = np.asarray(starts[chrom], dtype=np.int64)
        e = np.asarray(ends[chrom], dtype=np.int64)
        n = np.asarray(names[chrom], dtype=object)
        bt = np.asarray(biotypes[chrom], dtype=object)
        order = np.argsort(s, kind="stable")
        index[chrom] = {
            "start": s[order], "end": e[order],
            "name": n[order], "biotype": bt[order],
        }
    return index


def map_svs_to_genes(vcf_path: str, gene_index: dict, group: str):
    """Yield one dict per SV–gene overlap; returns (rows, n_sv, n_with_gene)."""
    rows = []
    n_sv = 0
    n_with_gene = 0

    vcf = cyvcf2.VCF(vcf_path)
    sample_idx = {s: vcf.samples.index(s) for s in SAMPLES if s in vcf.samples}

    for v in vcf:
        n_sv += 1
        pos = v.POS
        end_raw = v.INFO.get("END")
        end = int(end_raw) if end_raw is not None else pos
        sv_start = pos
        sv_end = max(end, pos)  # INS/BND: END == POS -> point overlap

        # Per-sample genotypes (positional over SAMPLES, "/" or "|" separated).
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
            "pos":     pos,
            "end":     end,
            "sv_id":   v.ID,
            "svtype":  v.INFO.get("SVTYPE"),
            "svlen":   v.INFO.get("SVLEN"),
            "support": v.INFO.get("SUPPORT"),
            "af":      v.INFO.get("AF"),
            **{f"gt_{s}": gts[s] for s in SAMPLES},
        }

        genes = gene_index.get(v.CHROM)
        if genes is None:
            continue

        mask = (genes["start"] <= sv_end) & (genes["end"] >= sv_start)
        hits = np.flatnonzero(mask)
        if hits.size == 0:
            continue
        n_with_gene += 1

        for k in hits:
            g_start = int(genes["start"][k])
            g_end = int(genes["end"][k])
            overlap_bp = min(sv_end, g_end) - max(sv_start, g_start) + 1
            rows.append({
                **base,
                "gene_id":      genes["name"][k],
                "gene_biotype": genes["biotype"][k],
                "gene_start":   g_start,
                "gene_end":     g_end,
                "overlap_bp":   overlap_bp,
            })

    return rows, n_sv, n_with_gene


def main():
    args = parse_args()

    group = os.path.basename(args.input).split(".")[0]
    os.makedirs(args.output, exist_ok=True)

    print(f"[sv_gene_mapping] Loading gene coordinates: {args.gff}", file=sys.stderr)
    gene_index = load_genes(args.gff)
    n_genes = sum(len(g["start"]) for g in gene_index.values())
    print(f"[sv_gene_mapping]   {n_genes:,} genes across {len(gene_index):,} contigs",
          file=sys.stderr)

    print(f"[sv_gene_mapping] Mapping SVs: {args.input} (group={group})", file=sys.stderr)
    rows, n_sv, n_with_gene = map_svs_to_genes(args.input, gene_index, group)

    df = pd.DataFrame(rows, columns=OUTPUT_COLS)
    out_tsv = os.path.join(args.output, f"{group}.sv_genes.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"[sv_gene_mapping]   SVs read           : {n_sv:,}", file=sys.stderr)
    print(f"[sv_gene_mapping]   SVs with gene hit  : {n_with_gene:,}", file=sys.stderr)
    print(f"[sv_gene_mapping]   SVs intergenic     : {n_sv - n_with_gene:,}", file=sys.stderr)
    print(f"[sv_gene_mapping]   SV–gene rows       : {len(df):,}", file=sys.stderr)
    print(f"[sv_gene_mapping] Written → {out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
