#!/usr/bin/env python3
"""
Build a parallel SV candidate table from the SnpEff-annotated multi-sample
Sniffles2 VCF.

This is a STANDALONE track: it does NOT feed the SNV+methylation genomic_score.
Each gene-associated SV becomes one row per affected gene, carrying per-sample
genotypes, a TAA-phenotype segregation pattern, and one-directional
cross-reference columns into the multi-omics ranked gene lists
(results/ranked_genes_lists/*.multiomics_ranked.tsv).

"Gene-associated" = the SV has a SnpEff impact of HIGH / MODERATE / LOW on the
gene (i.e. it actually overlaps the transcript/exon/structure). MODIFIER-only
associations (up/downstream/intergenic) are excluded to keep the table focused;
the full annotated SV VCF is retained on disk for the wider picture.

Reuses load_sorbi_to_loc() from build_ranked_genes.py so the gene_label key
(NCBI LOC* symbol) matches the ranked lists exactly.

    ./docker/run.sh python3 workflow/scripts/build_sv_table.py \
        --vcf        results/sv_calling/combined.annotated.vcf.gz \
        --gene-info  resources/NCBI_FTP/gene_info_4558 \
        --ranked-dir results/ranked_genes_lists \
        --out        results/sv_candidates/sv_candidate_table.tsv
"""

import argparse
import glob
import os
import sys

import cyvcf2
import pandas as pd

from build_ranked_genes import load_sorbi_to_loc

IMPACT_ORDER  = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
GENE_IMPACTS  = {"HIGH", "MODERATE", "LOW"}

# Phenotype groups (CLAUDE.md — TAA concentration in juice).
SAMPLES   = ["SBC4", "SBC10", "SBC11", "SBC23"]
TAA_HIGH  = ["SBC4", "SBC11", "SBC23"]
TAA_LOW   = ["SBC10"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",        required=True, help="SnpEff-annotated combined SV VCF (.vcf.gz)")
    p.add_argument("--gene-info",  required=True, help="NCBI gene_info tab file")
    p.add_argument("--ranked-dir", required=True, help="Directory of *.multiomics_ranked.tsv files")
    p.add_argument("--out",        required=True, help="Output TSV path")
    return p.parse_args()


def gt_str(g) -> str:
    """cyvcf2 genotype [a1, a2, phased] -> 0/0, 0|1, ./. etc."""
    a1, a2, phased = g[0], g[1], g[2]
    sep = "|" if phased else "/"
    tok = lambda a: "." if a is None or a < 0 else str(a)
    return f"{tok(a1)}{sep}{tok(a2)}"


def is_present(g) -> bool:
    """True if the genotype carries an alternate allele."""
    return (g[0] is not None and g[0] > 0) or (g[1] is not None and g[1] > 0)


def taa_pattern(present: dict) -> str:
    """Classify how an SV segregates with the TAA phenotype."""
    if not any(present.values()):
        return "none"
    high_all_present = all(present[s] for s in TAA_HIGH)
    high_all_absent  = all(not present[s] for s in TAA_HIGH)
    low_all_present  = all(present[s] for s in TAA_LOW)
    low_all_absent   = all(not present[s] for s in TAA_LOW)
    if high_all_present and low_all_absent:
        return "high-specific"
    if low_all_present and high_all_absent:
        return "low-specific"
    return "mixed"


def gene_worst_impacts(variant, sorbi_to_loc: dict) -> dict:
    """{gene_label: worst_impact} for genes this SV affects (HIGH/MODERATE/LOW)."""
    gene_worst: dict[str, str] = {}
    ann = variant.INFO.get("ANN")
    if not ann:
        return gene_worst
    for entry in str(ann).split(","):
        parts = entry.split("|")
        if len(parts) < 4:
            continue
        impact, sorbi_id = parts[2], parts[3]
        if impact not in GENE_IMPACTS:
            continue
        if not sorbi_id.startswith("SORBI_"):
            continue
        loc = sorbi_to_loc.get(sorbi_id)
        if not loc:
            continue
        prev = gene_worst.get(loc)
        if prev is None or IMPACT_ORDER.index(impact) < IMPACT_ORDER.index(prev):
            gene_worst[loc] = impact
    return gene_worst


def load_ranked(ranked_dir: str) -> dict:
    """gene_label -> list of (sample, rank, variant_score, methylation_score).

    Rank is the 1-based row position (files are pre-sorted by the ranking step).
    """
    ref: dict[str, list] = {}
    paths = sorted(glob.glob(os.path.join(ranked_dir, "*.multiomics_ranked.tsv")))
    for path in paths:
        sample = os.path.basename(path).split(".")[0]
        df = pd.read_csv(path, sep="\t")
        for rank, row in enumerate(df.itertuples(index=False), start=1):
            ref.setdefault(row.gene_label, []).append(
                (sample, rank, float(row.variant_score), float(row.methylation_score))
            )
    print(f"[sv_table] Ranked lists loaded: {len(paths)} ({[os.path.basename(p).split('.')[0] for p in paths]})",
          file=sys.stderr)
    return ref


def crossref(gene: str, ranked: dict) -> dict:
    hits = ranked.get(gene)
    if not hits:
        return {
            "in_ranked_list": "no",
            "ranked_in_samples": "",
            "best_variant_score": "",
            "best_methylation_score": "",
            "best_rank": "",
        }
    return {
        "in_ranked_list": "yes",
        "ranked_in_samples": ",".join(sorted({h[0] for h in hits})),
        "best_variant_score": round(max(h[2] for h in hits), 4),
        "best_methylation_score": round(max(h[3] for h in hits), 4),
        "best_rank": min(h[1] for h in hits),
    }


def main():
    args = parse_args()

    print(f"[sv_table] Loading gene ID map …", file=sys.stderr)
    sorbi_to_loc = load_sorbi_to_loc(args.gene_info)

    ranked = load_ranked(args.ranked_dir)

    print(f"[sv_table] Parsing SV VCF: {args.vcf}", file=sys.stderr)
    vcf = cyvcf2.VCF(args.vcf)
    vcf_samples = list(vcf.samples)
    sample_idx = {s: i for i, s in enumerate(vcf_samples)}

    rows = []
    n_sv = 0
    for variant in vcf:
        n_sv += 1
        genotypes = variant.genotypes  # list aligned with vcf_samples

        gts, present = {}, {}
        for s in SAMPLES:
            if s in sample_idx:
                g = genotypes[sample_idx[s]]
                gts[s] = gt_str(g)
                present[s] = is_present(g)
            else:
                gts[s] = "./."
                present[s] = False

        gene_worst = gene_worst_impacts(variant, sorbi_to_loc)
        if not gene_worst:
            continue  # not gene-associated → excluded from the candidate table

        pattern = taa_pattern(present)
        svtype  = variant.INFO.get("SVTYPE")
        svlen   = variant.INFO.get("SVLEN")
        end     = variant.INFO.get("END")
        support = variant.INFO.get("SUPPORT")
        end     = end if end is not None else variant.end

        base = {
            "chr": variant.CHROM,
            "pos": variant.POS,
            "end": end,
            "svtype": svtype if svtype is not None else "",
            "svlen": svlen if svlen is not None else "",
            "support": support if support is not None else "",
            "GT_SBC4": gts["SBC4"],
            "GT_SBC10": gts["SBC10"],
            "GT_SBC11": gts["SBC11"],
            "GT_SBC23": gts["SBC23"],
            "taa_pattern": pattern,
        }

        for gene, impact in gene_worst.items():
            row = dict(base)
            row["gene_label"] = gene
            row["impact"] = impact
            row.update(crossref(gene, ranked))
            rows.append(row)

    columns = [
        "chr", "pos", "end", "svtype", "svlen", "support",
        "GT_SBC4", "GT_SBC10", "GT_SBC11", "GT_SBC23",
        "taa_pattern", "gene_label", "impact",
        "in_ranked_list", "ranked_in_samples",
        "best_variant_score", "best_methylation_score", "best_rank",
    ]
    df = pd.DataFrame(rows, columns=columns)

    print(f"[sv_table] SVs parsed: {n_sv:,}  |  gene-associated rows: {len(df):,}",
          file=sys.stderr)

    if not df.empty:
        # Surface: in ranked list → worst impact → phenotype-segregating → score.
        df["_ranked"] = (df["in_ranked_list"] != "yes").astype(int)
        df["_impact"] = df["impact"].map(lambda x: IMPACT_ORDER.index(x))
        df["_seg"]    = (~df["taa_pattern"].isin(["high-specific", "low-specific"])).astype(int)
        df["_vscore"] = pd.to_numeric(df["best_variant_score"], errors="coerce").fillna(0.0)
        df = (
            df.sort_values(
                ["_ranked", "_impact", "_seg", "_vscore"],
                ascending=[True, True, True, False],
            )
            .drop(columns=["_ranked", "_impact", "_seg", "_vscore"])
            .reset_index(drop=True)
        )
        n_xref = (df["in_ranked_list"] == "yes").sum()
        n_seg  = df["taa_pattern"].isin(["high-specific", "low-specific"]).sum()
        print(f"[sv_table] Cross-referenced into ranked lists: {n_xref:,} rows", file=sys.stderr)
        print(f"[sv_table] Phenotype-segregating rows         : {n_seg:,}", file=sys.stderr)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[sv_table] Written → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
