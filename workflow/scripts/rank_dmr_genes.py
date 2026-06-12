#!/usr/bin/env python3
"""
Rank genes from DMR_annotated.tsv by diff.Methy (descending).
Outputs two tables:
  1. Hypermethylated in <accession>
  2. Hypomethylated in <accession> (= hypermethylated in the comparison partner)
Only rows with a non-missing gene_label are included.

./docker/run.sh python3 analysis/scripts/rank_dmr_genes.py \
    --output-hyper analysis/data/ranked_genes_lists/EPI-hyper_SBC10_genes.tsv \
    --output-hypo  analysis/data/ranked_genes_lists/EPI-hypo_SBC10_genes.tsv
"""

import argparse
import pandas as pd

DEFAULT_INPUT = "analysis/data/taa_DMR/DMR_annotated.tsv"


def rank_genes(df: pd.DataFrame, direction_value: str, feature: str | None = "promoter") -> pd.DataFrame:
    mask = (df["direction"] == direction_value) & (df["gene_label"].notna()) & (df["gene_label"].str.strip() != "")
    if feature:
        mask &= df["feature"] == feature
    filtered = df[mask].copy()

    return (
        filtered
        .sort_values("diff.Methy", ascending=False)
        .drop_duplicates(subset="gene_label")
        [["gene_label", "diff.Methy", "direction", "feature", "chr", "start", "end"]]
        .reset_index(drop=True)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Path to DMR_annotated.tsv")
    parser.add_argument("--accession", default="SBC10",
                        help="Focal accession (default: SBC10)")
    parser.add_argument("--output-hyper", default=None,
                        help="Save hypermethylated table to this TSV path")
    parser.add_argument("--output-hypo", default=None,
                        help="Save hypomethylated table to this TSV path")
    parser.add_argument("--feature", default="promoter",
                        help="Filter to this genomic feature (default: promoter). Pass '' to include all.")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # Rows involving the focal accession
    acc = args.accession
    involved = df[(df["sample_a"] == acc) | (df["sample_b"] == acc)]

    hyper_direction = f"hyper_{acc}"
    # Hypomethylated in acc = the other sample is hyper in a comparison involving acc
    hypo_directions = involved.loc[involved["direction"] != hyper_direction, "direction"].unique()

    feature = args.feature or None
    hyper_table = rank_genes(involved, hyper_direction, feature)
    # Combine all hypo directions, then re-rank across them
    hypo_parts = [rank_genes(involved, d, feature) for d in hypo_directions]
    hypo_table = (
        pd.concat(hypo_parts)
        .sort_values("diff.Methy", ascending=True)
        .drop_duplicates(subset="gene_label")
        .reset_index(drop=True)
    )

    pd.set_option("display.max_rows", None)
    pd.set_option("display.width", 120)

    feature_label = f" [{feature}]" if feature else ""
    print(f"=== Hypermethylated in {acc}{feature_label} ({len(hyper_table)} genes) ===")
    if args.output_hyper:
        hyper_table.to_csv(args.output_hyper, sep="\t", index=False)
        print(f"Saved to {args.output_hyper}")
    else:
        print(hyper_table.to_string(index=False))

    print()

    print(f"=== Hypomethylated in {acc}{feature_label} ({len(hypo_table)} genes) ===")
    if args.output_hypo:
        hypo_table.to_csv(args.output_hypo, sep="\t", index=False)
        print(f"Saved to {args.output_hypo}")
    else:
        print(hypo_table.to_string(index=False))


if __name__ == "__main__":
    main()
