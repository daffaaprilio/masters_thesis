#!/usr/bin/env python3
"""Export the candidate gene lists that are too long to print in the appendix.

For each phenotype notebook (taa, dryness, more_juiciness) this writes two CSVs:

  * variant candidates  -- gene-level score >= 0.5 (the `genomics_cand` list the
    fixed-background ORA in each notebook is run on)
  * promoter-DMR candidates -- the `phenotype_candidates` table (promoter
    hypomethylated in the contrast accession in >= 1 phenotype contrast)

The scoring/parsing helpers are exec'd straight out of dryness.ipynb so this
script and the notebooks cannot drift apart.
"""

import json
from pathlib import Path

import cyvcf2
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analysis"
OUT_DIR = ROOT / "writing" / "supplementary"
DMR_PATH = ROOT / "results" / "DMR" / "DMR_annotated.tsv"
HELPER_NB, HELPER_CELL = ANALYSIS / "dryness.ipynb", 5

# phenotype -> (vcf, contrast accession, the accessions it is contrasted against)
PHENOTYPES = {
    "taa":             dict(vcf="SBC10.vcf.gz",     focal="SBC10", others=["SBC4", "SBC11", "SBC23"]),
    "dryness":         dict(vcf="SBC11.vcf.gz",     focal="SBC11", others=["SBC4", "SBC10", "SBC23"]),
    "extra_juiciness": dict(vcf="notSBC11.vcf.gz",  focal="SBC11", others=["SBC4", "SBC10", "SBC23"]),
}

VARIANT_COLS = ["ann_gene_id", "ann_gene_name", "score", "percentile", "rank",
                "ann_impact", "ann_effect", "ann_hgvs_p", "chrom", "pos",
                "sift_score", "sift_prediction"]

# CSV header names, chosen so both file types share `gene_id` and drop the
# internal `ann_`/`gene_label` prefixes used in the notebooks.
VARIANT_RENAME = {"ann_gene_id": "gene_id", "ann_gene_name": "gene_name",
                  "ann_impact": "impact", "ann_effect": "effect",
                  "ann_hgvs_p": "hgvs_p"}


def load_helpers():
    """exec the notebook cell that defines parse_vcf / scoring / merge_to_gene_max."""
    cell = json.load(open(HELPER_NB))["cells"][HELPER_CELL]
    ns = {"pd": pd, "np": np, "cyvcf2": cyvcf2, "Path": Path}
    exec("".join(cell["source"]), ns)
    return ns


def variant_candidates(ns, vcf_name):
    df = ns["parse_vcf"](str(ANALYSIS / "data" / vcf_name))
    # notebook cell: drop multi-gene ANN gene ids ("-" / "&" joined)
    df = df[(~df["ann_gene_id"].str.contains("-", regex=False, na=False, case=False))
            & (~df["ann_gene_id"].str.contains("&", regex=False, na=False, case=False))]
    df_genes = ns["merge_to_gene_max"](ns["scoring"](df))
    cand = df_genes[df_genes["score"] >= 0.5].copy()
    cols = [c for c in VARIANT_COLS if c in cand.columns]
    return (cand[cols].sort_values("score", ascending=False)
            .reset_index(drop=True).rename(columns=VARIANT_RENAME))


def dmr_candidates(df_dmr, focal, others):
    promoter = df_dmr[(df_dmr["feature"] == "promoter")
                      & df_dmr["gene_label"].notna()
                      & (df_dmr["gene_label"].astype(str).str.strip() != "")].copy()
    contrast = promoter[
        ((promoter["sample_a"] == focal) & promoter["sample_b"].isin(others))
        | ((promoter["sample_b"] == focal) & promoter["sample_a"].isin(others))
    ].copy()
    contrast["hypo_in_focal"] = contrast["direction"] != f"hyper_{focal}"
    return (contrast[contrast["hypo_in_focal"]]
            .groupby("gene_label")
            .agg(n_supporting_contrasts=("hypo_in_focal", "sum"),
                 mean_methyl_diff=("diff.Methy", "mean"))
            .sort_values("mean_methyl_diff", key=lambda s: s.abs(), ascending=False)
            .reset_index()
            .rename(columns={"gene_label": "gene_id"}))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ns = load_helpers()
    df_dmr = pd.read_csv(DMR_PATH, sep="\t")

    n = 0
    summary = []
    for pheno, cfg in PHENOTYPES.items():
        var = variant_candidates(ns, cfg["vcf"])
        n += 1
        p = OUT_DIR / f"S{n}_{pheno}_variant_candidates.csv"
        var.to_csv(p, index=False)
        summary.append((p.name, f"{pheno} / variant (score >= 0.5)", len(var)))

        dmr = dmr_candidates(df_dmr, cfg["focal"], cfg["others"])
        n += 1
        p = OUT_DIR / f"S{n}_{pheno}_dmr_candidates.csv"
        dmr.to_csv(p, index=False)
        summary.append((p.name, f"{pheno} / promoter DMR", len(dmr)))

    print(f"\n{'file':52s} {'contents':40s} {'genes':>7s}")
    print("-" * 102)
    for name, what, k in summary:
        print(f"{name:52s} {what:40s} {k:7,d}")


if __name__ == "__main__":
    main()
