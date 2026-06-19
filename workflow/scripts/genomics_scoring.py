#!/opt/conda/bin/python3

import os
import sys
import numpy as np
import cyvcf2
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")          # headless: write figures to file, never to a display
import matplotlib.pyplot as plt
import argparse
import gzip, re

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Reads VCF file, assigns scores to each gene containing variants listed in that VCF file.
        Variants are associated with the samples' phenotype, i.e., SBC10-private variants are associated with low TAA production (aconitate isomerase), SBC11-private variants; D gene (stem juiciness).
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input",
        help="Path to the input sample-exclusive VCF file",
    )

    parser.add_argument(
        "-o", "--output", 
        help="Parent directory of the output files",
    )

    parser.add_argument(
        "-g", "--gff",
        help="Path to the GFF3 file for gene length correction reference",
    )

    return parser.parse_args()

# prepare dataframe
ANN_FIELDS = [
    "ann_allele", "ann_effect", "ann_impact", "ann_gene_name", "ann_gene_id",
    "ann_feature_type", "ann_feature_id", "ann_biotype", "ann_rank",
    "ann_hgvs_c", "ann_hgvs_p", "ann_cdna_pos", "ann_cds_pos", "ann_aa_pos",
    "ann_distance", "ann_extra",
]

SIFT_FIELDS = [
    "sift_allele", "sift_transcript", "sift_gene_id", "sift_gene_name",
    "sift_region", "sift_variant_type", "sift_aa_change", "sift_aa_pos",
    "sift_score", "sift_median", "sift_num_seqs", "sift_allele_type",
    "sift_prediction",
]

# SIFT4G writes DELETERIOUS (score < 0.05) or TOLERATED; NA means non-coding
_SIFT_PRIORITY = {"DELETERIOUS": 0, "TOLERATED": 1}

# impact severity for tiebreaking (HIGH most severe -> lowest rank)
_IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def _parse_ann(raw: str) -> list[dict]:
    records = []
    for entry in raw.split(","):
        parts = entry.split("|")
        parts += [""] * (len(ANN_FIELDS) - len(parts))
        records.append(dict(zip(ANN_FIELDS, parts[:len(ANN_FIELDS)])))
    return records


def _parse_siftinfo(raw: str) -> dict:
    """Returns (allele, transcript) -> best SIFT record.

    Both SnpEff and SIFT4G are built from NCBIv3 GFF, so ann_feature_id and
    sift_transcript share the same XM_* namespace and can be joined directly.
    SIFT4G annotates one representative transcript per gene; unmatched ANN
    transcript rows receive no SIFT columns.
    """
    by_key: dict[tuple, dict] = {}
    for entry in raw.split(","):
        parts = entry.split("|")
        parts += [""] * (len(SIFT_FIELDS) - len(parts))
        d = dict(zip(SIFT_FIELDS, parts[:len(SIFT_FIELDS)]))
        key = (d["sift_allele"], d["sift_transcript"])
        pred = d.get("sift_prediction", "")
        current = by_key.get(key)
        if current is None or _SIFT_PRIORITY.get(pred, 99) < _SIFT_PRIORITY.get(current.get("sift_prediction", ""), 99):
            by_key[key] = d
    return by_key


def _reduce_annotations(annotations: list[dict], sift_by_key: dict) -> list[tuple]:
    """Reduce a variant's many ANN entries to ONE per gene (see MULT_ANN_TO_ONE_ANN.md).

    SnpEff emits one ANN per affected transcript, but SIFT4G anchors to a single
    representative transcript per gene. To make the two comparable we keep one
    ANN per gene, chosen by a single deterministic sort (lowest key wins):

      1. SIFT match  -- (allele, transcript) is a SIFTINFO key. This is Case 1
         (apple-to-apple): it pairs the SnpEff effect and the SIFT score on the
         SAME transcript. Top key, so a SIFT-anchored row always wins when present.

    When no entry matches SIFT (Case 2, gene-level umbrella), keys 2-4 decide:
      2. impact severity  -- HIGH > MODERATE > LOW > MODIFIER
      3. SnpEff ANN order -- severity proxy within an impact tier, captured here
         while the original ANN order is still intact
      4. biotype          -- prefer protein_coding over pseudogene / nc

    Returns a list of (ann, sift) pairs, one per gene.
    """
    by_gene: dict[str, list[tuple]] = {}
    for order, ann in enumerate(annotations):
        key = (ann.get("ann_allele", ""), ann.get("ann_feature_id", ""))
        sort_key = (
            0 if key in sift_by_key else 1,                          # 1. SIFT match (Case 1)
            _IMPACT_RANK.get(ann.get("ann_impact", ""), 99),         # 2. impact severity
            order,                                                   # 3. SnpEff ANN order
            0 if ann.get("ann_biotype") == "protein_coding" else 1,  # 4. biotype
        )
        by_gene.setdefault(ann.get("ann_gene_id", ""), []).append((sort_key, key, ann))

    reduced = []
    for entries in by_gene.values():
        _, key, ann = min(entries, key=lambda t: t[0])              # winner of the sort
        reduced.append((ann, sift_by_key.get(key, {})))
    return reduced


def parse_vcf(vcf_path: str) -> pd.DataFrame:
    """
    Parse a SIFT4G+SnpEff annotated VCF into a DataFrame.

    One row per (variant, gene): the many SnpEff ANN transcript entries are
    reduced to a single representative transcript per gene via
    _reduce_annotations (see MULT_ANN_TO_ONE_ANN.md), and the SIFT columns of
    that transcript are joined on. This is the key difference from the baseline
    notebook, which keeps one row per ANN transcript entry.
    """
    vcf = cyvcf2.VCF(vcf_path)
    rows = []
    for variant in vcf:
        gt_arr = variant.genotypes[0]
        sep = "|" if gt_arr[2] else "/"
        base = {
            "chrom":  variant.CHROM,
            "pos":    variant.POS,
            "ref":    variant.REF,
            "alt":    ",".join(variant.ALT),
            "qual":   variant.QUAL,
            "filter": variant.FILTER,
            "gt":     f"{gt_arr[0]}{sep}{gt_arr[1]}",
        }

        ann_raw  = variant.INFO.get("ANN")
        sift_raw = variant.INFO.get("SIFTINFO")
        if not ann_raw:
            rows.append(base)
            continue

        annotations = _parse_ann(ann_raw)
        sift_by_key = _parse_siftinfo(sift_raw) if sift_raw else {}
        for ann, sift in _reduce_annotations(annotations, sift_by_key):
            rows.append({**base, **ann, **sift})

    df = pd.DataFrame(rows)
    df["sift_score"] = pd.to_numeric(df.get("sift_score"), errors="coerce")
    df["sift_score_c"] = 1 - df["sift_score"]
    return df

def arbitrary_score(df):
    """
    Obtain arbitrary score for rows with SnpEff input but no SIFT output
    For HIGH, though, ignore the average SIFT score.
    SIFT score for HIGH SnpEff sites are misleadingly low, due to NA-by-0-imputation
    """
    scores = df["sift_score_c"].fillna(0)
    return scores.groupby(df["ann_impact"]).mean().to_dict()

def scoring(df):
    impact_default = arbitrary_score(df)
    score = df["ann_impact"].map(impact_default)
    score = score.mask(df["sift_score_c"].notna(), df["sift_score_c"])
    score = score.mask(df["ann_impact"] == "HIGH", 1.0)
    df["score"] = score
    return df

def merge_to_gene_max(df_scored: pd.DataFrame, gene_key: str = "ann_gene_id") -> pd.DataFrame:
    """Collapse variant/transcript rows to one row per gene.

    Maximum / worst-variant approach: each gene keeps the single annotation
    with the highest ``score``, so the gene-level risk equals its most
    damaging variant. The full worst-variant row is retained (effect, impact,
    SIFT, position) for context. Rows with no gene id are dropped.
    """
    scored = df_scored[df_scored[gene_key].fillna("").ne("")]
    idx = scored.groupby(gene_key)["score"].idxmax()       # worst variant per gene
    return (
        scored.loc[idx]
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )

def merge_to_gene_sum(df_scored: pd.DataFrame, gene_key: str = "ann_gene_id") -> pd.DataFrame:
    """Collapse variant/transcript rows to one row per gene by summing scores.

    Sum approach: a gene's risk is the total of every annotation score it
    carries, so it rewards genes hit by many and/or more damaging variants
    (unlike the max approach, which only looks at the single worst variant).
    n_rows is kept alongside for context. Rows with no gene id are dropped.

    Note: SnpEff emits one row per transcript, so multi-transcript genes
    accumulate more rows; the sum reflects that transcript multiplicity.
    """
    scored = df_scored[df_scored[gene_key].fillna("").ne("")]
    return (
        scored.groupby(gene_key)
        .agg(
            score=("score", "sum"),
            n_rows=("score", "size"),
            ann_gene_name=("ann_gene_name", "first"),
        )
        .reset_index()
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )

def load_gene_length(gff_path: str) -> dict:
    """Map gene Name (LOC...) -> genomic length in bp, from an NCBI GFF3.

    Parses 'gene' and 'pseudogene' features; length = end - start + 1 (genomic
    span, introns + UTRs). Keyed by the Name= attribute, which equals
    ann_gene_id for LOC-named genes. tRNA (SnpEff '..._N') and combo/intergenic
    keys won't match -> those genes get NaN length downstream.
    """
    lengths = {}
    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = re.search(r"Name=([^;]+)", f[8])
            if m:
                lengths[m.group(1)] = int(f[4]) - int(f[3]) + 1
    return lengths

def merge_to_gene_max_norm(df_scored: pd.DataFrame, gene_length: dict,
                           gene_key: str = "ann_gene_id", per_kb: bool = True) -> pd.DataFrame:
    """MAX gene score divided by genomic gene length (default: score per kb).

    Reuses merge_to_gene_max (full worst-variant row), then normalizes. The
    normalized value is named ``score`` (so plot_kdeplot works); the raw max is
    kept as score_raw. Genes with no length match (tRNA '_N', combo keys) ->
    NaN score.

    Note: max / length mostly re-ranks equal-max genes inversely by length.
    """
    g = merge_to_gene_max(df_scored, gene_key).rename(columns={"score": "score_raw"})
    g["gene_length"] = g[gene_key].map(gene_length)
    g["score"] = g["score_raw"] / g["gene_length"] * (1000 if per_kb else 1)
    return g.sort_values("score", ascending=False).reset_index(drop=True)

def merge_to_gene_sum_norm(df_scored: pd.DataFrame, gene_length: dict,
                           gene_key: str = "ann_gene_id", per_kb: bool = True) -> pd.DataFrame:
    """SUM gene score divided by genomic gene length (default: score per kb).

    Reuses merge_to_gene_sum, then normalizes. The normalized value is named
    ``score`` (so plot_kdeplot works); the raw sum is kept as score_raw.
    Genes with no length match (tRNA '_N', combo keys) -> NaN score.

    Note: length-normalization corrects for gene size only, NOT the SUM
    method's transcript-multiplicity inflation (see n_rows).
    """
    g = merge_to_gene_sum(df_scored, gene_key).rename(columns={"score": "score_raw"})
    g["gene_length"] = g[gene_key].map(gene_length)
    g["score"] = g["score_raw"] / g["gene_length"] * (1000 if per_kb else 1)
    return g.sort_values("score", ascending=False).reset_index(drop=True)

def plot_kdeplot(df: pd.DataFrame, title: str, ax=None, clip=None):
    """KDE of the per-row/per-gene ``score`` column (mirrors the notebook helper).

    The x-axis starts at 0.129 (the LOW impact floor) so the dense pile-up of
    near-zero MODIFIER scores doesn't flatten the rest of the distribution.
    """
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(8, 5))

    x_min = 0.129
    x_max = clip[1] if clip is not None else df["score"].max()
    effective_clip = (x_min, x_max)

    sns.kdeplot(
        data=df,
        x="score",
        fill=True,
        clip=effective_clip,
        cut=0,
        bw_adjust=1.0,
        ax=ax,
    )

    ax.set_xlim(effective_clip)
    ax.set_title(f"Scores distribution on {title} level")
    ax.set_xlabel("Score")
    ax.set_ylabel("Density")

    if own_fig:
        fig.tight_layout()
    return ax


def plotting_summary(gene_tables: dict, sample: str, out_path: str) -> str:
    """Render the 2x2 KDE comparison of the four gene-level scoring approaches.

    ``gene_tables`` maps {"max", "sum", "max_norm", "sum_norm"} -> per-gene
    DataFrame (each with a ``score`` column). Each panel is auto-scaled to its
    own range. Saves the figure to ``out_path`` and returns that path.
    """
    panels = [
        (gene_tables["max"],      "gene-MAX",      None),
        (gene_tables["sum"],      "gene-SUM",      None),
        (gene_tables["max_norm"], "gene-MAX / kb", None),
        (gene_tables["sum_norm"], "gene-SUM / kb", None),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for (d, t, c), ax in zip(panels, axes.flat):
        plot_kdeplot(d, t, ax=ax, clip=c)
    fig.suptitle(f"{sample}: gene-level score distributions", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    args = parse_args()
    sns.set_theme(style="whitegrid")

    sample = os.path.basename(args.input).split(".")[0]
    os.makedirs(args.output, exist_ok=True)

    print(f"[genomics_scoring] Parsing VCF: {args.input}", file=sys.stderr)
    df = parse_vcf(args.input)
    print(f"[genomics_scoring]   {len(df):,} (variant, gene) rows", file=sys.stderr)

    df_scored = scoring(df)

    print(f"[genomics_scoring] Loading gene lengths: {args.gff}", file=sys.stderr)
    gene_length = load_gene_length(args.gff)
    print(f"[genomics_scoring]   {len(gene_length):,} genes with length", file=sys.stderr)

    gene_tables = {
        "max":      merge_to_gene_max(df_scored),
        "sum":      merge_to_gene_sum(df_scored),
        "max_norm": merge_to_gene_max_norm(df_scored, gene_length),
        "sum_norm": merge_to_gene_sum_norm(df_scored, gene_length),
    }

    for name, table in gene_tables.items():
        out_tsv = os.path.join(args.output, f"{sample}.gene_{name}.tsv")
        table.to_csv(out_tsv, sep="\t", index=False)
        print(f"[genomics_scoring] Written → {out_tsv}  ({len(table):,} genes)", file=sys.stderr)

    fig_path = os.path.join(args.output, f"{sample}.scores_summary.png")
    plotting_summary(gene_tables, sample, fig_path)
    print(f"[genomics_scoring] Written → {fig_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
