#!/usr/bin/env python3
"""
Build a per-sample multi-omics ranked gene list.

Parses the SIFT-annotated VCF (which carries both SnpEff ANN and SIFTINFO
tags) to derive per-gene genomic disruption scores, then merges with the
annotated DMR table to obtain per-gene epigenomic scores.

Genomic score  = impact_score  +  0.5 * sift_disruption
  impact_score:    HIGH=1.0  MODERATE=0.67  LOW=0.33  MODIFIER=0.0
  sift_disruption: 1 - min(SIFT_score) across all scored CDS variants
                   (0 when no SIFT-scored variants exist for that gene)

Epigenomic score = max(|diff.Methy|) across all DMR rows involving the sample.

Output is sorted by genomic_score DESC, epigenomic_score DESC.
"""

import argparse
import sys

import cyvcf2
import pandas as pd

IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
IMPACT_SCORE = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",       required=True, help="SIFT-annotated VCF (.vcf.gz)")
    p.add_argument("--dmr",       required=True, help="DMR_annotated.tsv")
    p.add_argument("--gene-info", required=True, help="NCBI gene_info tab file")
    p.add_argument("--sample",    required=True, help="Sample name (e.g. SBC10)")
    p.add_argument("--out",       required=True, help="Output TSV path")
    return p.parse_args()


def load_sorbi_to_loc(gene_info_path: str) -> dict:
    gi = pd.read_csv(gene_info_path, sep="\t", usecols=["Symbol", "LocusTag"])
    return gi.set_index("LocusTag")["Symbol"].to_dict()


def parse_vcf(vcf_path: str, sorbi_to_loc: dict):
    """
    Returns:
        snpeff_hits : {loc_id: [impact, ...]}   — all SnpEff impacts seen for gene
        sift_hits   : {loc_id: [sift_score, ...]} — SIFT scores for scored CDS variants
    """
    snpeff_hits: dict[str, list] = {}
    sift_hits:   dict[str, list] = {}

    vcf = cyvcf2.VCF(vcf_path)
    for variant in vcf:

        # ── SnpEff ANN field ───────────────────────────────────────────────
        ann = variant.INFO.get("ANN")
        if ann:
            for entry in str(ann).split(","):
                parts = entry.split("|")
                if len(parts) < 4:
                    continue
                impact   = parts[2]
                sorbi_id = parts[3]
                if impact not in IMPACT_SCORE:
                    continue
                if not sorbi_id.startswith("SORBI_"):
                    continue
                loc_id = sorbi_to_loc.get(sorbi_id)
                if loc_id:
                    snpeff_hits.setdefault(loc_id, []).append(impact)

        # ── SIFT4G SIFTINFO field ──────────────────────────────────────────
        # Format: Allele|Transcript|GeneId|GeneName|Region|VariantType|
        #         RefAA/AltAA|AminoPos|SIFT_score|SIFT_median|NUM_seqs|
        #         Allele_Type|SIFT_prediction
        siftinfo = variant.INFO.get("SIFTINFO")
        if siftinfo:
            for entry in str(siftinfo).split(","):
                parts = entry.split("|")
                if len(parts) < 13:
                    continue
                loc_id    = parts[2]
                sift_pred = parts[12]
                score_str = parts[8]
                if sift_pred not in ("TOLERATED", "DAMAGING"):
                    continue
                try:
                    sift_hits.setdefault(loc_id, []).append(float(score_str))
                except ValueError:
                    pass

    return snpeff_hits, sift_hits


def build_genomics(snpeff_hits: dict, sift_hits: dict) -> pd.DataFrame:
    all_genes = set(snpeff_hits) | set(sift_hits)
    rows = []
    for gene in all_genes:
        impacts = snpeff_hits.get(gene, [])
        scores  = sift_hits.get(gene, [])

        worst_impact = (
            min(impacts, key=lambda x: IMPACT_ORDER.index(x))
            if impacts else "MODIFIER"
        )
        impact_score = IMPACT_SCORE[worst_impact]

        if scores:
            min_sift      = min(scores)
            sift_disruption = 1.0 - min_sift
        else:
            min_sift        = None
            sift_disruption = 0.0

        rows.append({
            "gene_label":      gene,
            "worst_impact":    worst_impact,
            "impact_score":    round(impact_score, 4),
            "min_sift_score":  round(min_sift, 4) if min_sift is not None else None,
            "sift_disruption": round(sift_disruption, 4),
            "genomic_score":   round(impact_score + sift_disruption, 4),
        })

    return pd.DataFrame(rows)


def build_epigenomics(dmr_path: str, sample: str) -> pd.DataFrame:
    dmr = pd.read_csv(dmr_path, sep="\t")

    # Keep rows where the sample participates in either column
    dmr = dmr[
        (dmr["sample_a"] == sample) | (dmr["sample_b"] == sample)
    ].copy()
    dmr = dmr.dropna(subset=["gene_label"])
    dmr = dmr[dmr["gene_label"].astype(str).str.strip() != ""]
    dmr["abs_methy"] = dmr["diff.Methy"].abs()

    # Per-gene: take the DMR row with the strongest methylation difference
    idx_max = dmr.groupby("gene_label")["abs_methy"].idxmax()
    epi = dmr.loc[idx_max, ["gene_label", "abs_methy", "direction", "feature"]].copy()
    epi = epi.rename(columns={"abs_methy": "epigenomic_score"})
    epi["epigenomic_score"] = epi["epigenomic_score"].round(4)
    return epi.reset_index(drop=True)


def main():
    args = parse_args()

    print(f"[ranked_genes] Loading gene ID map …", file=sys.stderr)
    sorbi_to_loc = load_sorbi_to_loc(args.gene_info)

    print(f"[ranked_genes] Parsing VCF: {args.vcf}", file=sys.stderr)
    snpeff_hits, sift_hits = parse_vcf(args.vcf, sorbi_to_loc)
    print(
        f"[ranked_genes]   SnpEff-mapped genes : {len(snpeff_hits):,}",
        f"\n[ranked_genes]   SIFT-scored genes   : {len(sift_hits):,}",
        file=sys.stderr,
    )

    genomics = build_genomics(snpeff_hits, sift_hits)
    print(f"[ranked_genes] Genomics table: {len(genomics):,} genes", file=sys.stderr)

    print(f"[ranked_genes] Parsing DMR: {args.dmr} (sample={args.sample})", file=sys.stderr)
    epigenomics = build_epigenomics(args.dmr, args.sample)
    print(f"[ranked_genes] Epigenomics table: {len(epigenomics):,} genes", file=sys.stderr)

    merged = genomics.merge(epigenomics, on="gene_label", how="outer")
    merged["genomic_score"]    = merged["genomic_score"].fillna(0.0)
    merged["epigenomic_score"] = merged["epigenomic_score"].fillna(0.0)
    merged["worst_impact"]     = merged["worst_impact"].fillna("MODIFIER")
    merged["impact_score"]     = merged["impact_score"].fillna(0.0)
    merged["sift_disruption"]  = merged["sift_disruption"].fillna(0.0)

    both = (merged["genomic_score"] > 0) & (merged["epigenomic_score"] > 0)
    print(
        f"[ranked_genes] Genes in both layers : {both.sum():,}",
        f"\n[ranked_genes] Total genes          : {len(merged):,}",
        file=sys.stderr,
    )

    max_g = merged["genomic_score"].max()
    max_e = merged["epigenomic_score"].max()
    merged["variant_score"]     = (merged["genomic_score"]    / max_g).round(4)
    merged["methylation_score"] = (merged["epigenomic_score"] / max_e).round(4)

    merged = merged.sort_values(
        ["variant_score", "methylation_score"], ascending=False
    ).reset_index(drop=True)

    merged[["gene_label", "variant_score", "methylation_score"]].to_csv(
        args.out, sep="\t", index=False
    )
    print(f"[ranked_genes] Written → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
