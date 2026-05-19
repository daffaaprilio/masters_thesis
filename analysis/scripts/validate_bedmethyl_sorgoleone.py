#!/usr/bin/env python3
"""
Validation script for filtered bedMethyl files from sorgoleone loci.

Runs 7 QC checks and writes a summary TSV to analysis/qc/.

Usage:
    python analysis/scripts/validate_sorgoleone_bedmethyl.py \
        [--indir analysis/data/sorgoleone_bedmethyl] \
        [--outdir analysis/qc] \
        [--gff resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff]
"""

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

BEDMETHYL_COLS = [
    "chrom", "start", "end", "code", "coverage", "strand",
    "start2", "end2", "color", "valid_coverage", "frac_modified",
    "n_mod", "n_canonical", "n_other_mod", "n_delete", "n_fail", "n_diff", "n_nocall",
]
GFF3_COLS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

ARTEFACT_CODE = "21839"
PRIMARY_CODE = "m"    # 5mC
SECONDARY_CODE = "a"  # 6mA

LOC_LABELS = {
    "LOC8066368":   "SbDES2",
    "LOC8079957":   "SbDES3",
    "LOC8080259":   "SbOMT3",
    "LOC8081692":   "SbCYP71AM1",
    "LOC110435045": "SbDES2-homolog",
    "LOC8072903":   "SbDES3-homolog-1",
    "LOC8055482":   "SbDES3-homolog-2",
    "LOC8079958":   "SbDES3-homolog-3",
    "LOC8076922":   "SbOMT3-homolog-1",
    "LOC110436225": "SbOMT3-homolog-2",
    "LOC8085153":   "SbOMT3-homolog-3",
}

WARNINGS = []


def setup_logging(log_path):
    """Configure root logger to write to console (clean) and log file (timestamped)."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter("%(message)s"))
    fh = logging.FileHandler(log_path)
    fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                                      datefmt="%Y-%m-%d %H:%M:%S"))
    root.addHandler(ch)
    root.addHandler(fh)


def flag(msg):
    WARNINGS.append(msg)
    logging.warning(f"*** {msg}")


def section(title):
    bar = "=" * 62
    logging.info(f"\n{bar}\n  {title}\n{bar}")


def parse_gff_attribute(attrs, key):
    for field in attrs.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def load_gene_regions(gff_path, flank):
    gff = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=GFF3_COLS)
    genes = gff[gff["feature"] == "gene"].copy()
    genes["loc_name"] = genes["attributes"].apply(lambda a: parse_gff_attribute(a, "Name"))
    sel = genes[genes["loc_name"].isin(LOC_LABELS)].copy()
    sel["label"] = sel["loc_name"].map(LOC_LABELS)
    sel["region_start"] = (sel["start"] - 1 - flank).clip(lower=0).astype(int)
    sel["region_end"] = (sel["end"] + flank).astype(int)
    return sel[["seqname", "region_start", "region_end", "label"]].reset_index(drop=True)


def assign_locus(df, regions):
    """Assign each site to a gene locus; unassigned if outside all regions."""
    result = pd.Series("unassigned", index=df.index, dtype=str)
    for _, reg in regions.iterrows():
        mask = (
            (df["chrom"] == reg["seqname"]) &
            (df["start"] >= reg["region_start"]) &
            (df["end"] <= reg["region_end"])
        )
        result[mask] = reg["label"]
    return result


def load_sample(path):
    return pd.read_csv(
        path, sep="\t", header=None, names=BEDMETHYL_COLS,
        dtype={"code": str, "valid_coverage": float, "frac_modified": float},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indir", default="analysis/data/sorgoleone_bedmethyl",
                        help="Directory containing .sorgoleone.bed files")
    parser.add_argument("--outdir", default="analysis/data/sorgoleone_bedmethyl",
                        help="Output directory for the validation TSV and log")
    parser.add_argument("--gff", default=None,
                        help="GFF3 annotation file for per-locus site counts (Check 2)")
    parser.add_argument("--flank", type=int, default=2000,
                        help="Flank used during filtering, for GFF region extraction (default: 2000)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = Path(__file__).parent.parent / "logs" / f"validate_bedmethyl_sorgoleone_{timestamp}.log"
    setup_logging(log_path)

    logging.info(f"Log: {log_path}")
    logging.info(f"Input dir: {args.indir}")

    indir = Path(args.indir)
    bed_files = sorted(indir.glob("*.sorgoleone.bed"))
    if not bed_files:
        logging.error(f"No .sorgoleone.bed files found in {indir}")
        sys.exit(1)

    samples = {p.name.split(".")[0]: p for p in bed_files}
    logging.info(f"Samples: {', '.join(samples)}")

    raw  = {s: load_sample(p) for s, p in samples.items()}
    data = {s: df[df["code"] != ARTEFACT_CODE].copy() for s, df in raw.items()}
    mc   = {s: df[df["code"] == PRIMARY_CODE].copy() for s, df in data.items()}

    tsv_rows = []

    # ------------------------------------------------------------------ #
    # Check 1: Per-file site counts
    # ------------------------------------------------------------------ #
    section("CHECK 1: Per-file site counts")
    for s, df in data.items():
        n_m = (df["code"] == PRIMARY_CODE).sum()
        n_a = (df["code"] == SECONDARY_CODE).sum()
        total = len(df)
        logging.info(f"  {s}: total={total:,}  5mC={n_m:,}  6mA={n_a:,}")
        if total == 0:
            flag(f"{s}: zero sites after artefact exclusion")
        tsv_rows += [
            {"sample": s, "check": "site_count_total", "value": total},
            {"sample": s, "check": "site_count_5mC",   "value": int(n_m)},
            {"sample": s, "check": "site_count_6mA",   "value": int(n_a)},
        ]

    # ------------------------------------------------------------------ #
    # Check 2: Per-locus 5mC site counts
    # ------------------------------------------------------------------ #
    section("CHECK 2: Per-locus 5mC site counts")
    if args.gff:
        regions = load_gene_regions(args.gff, args.flank)
        logging.info(f"  Loaded {len(regions)}/{len(LOC_LABELS)} gene regions from GFF3.")
        locus_counts = {}
        for s, df in mc.items():
            df = df.copy()
            df["locus"] = assign_locus(df, regions)
            locus_counts[s] = df.groupby("locus").size()

        col_w = max(len(s) for s in samples)
        header = f"  {'Locus':<26}" + "".join(f"{s:>{col_w + 2}}" for s in samples)
        logging.info(header)
        for label in sorted(LOC_LABELS.values()):
            row = f"  {label:<26}"
            for s in samples:
                n = int(locus_counts[s].get(label, 0))
                row += f"{n:>{col_w + 2},}"
                tsv_rows.append({"sample": s, "check": f"locus_5mC_{label}", "value": n})
                if n == 0:
                    flag(f"{s}: 0 5mC sites at {label} — possible filtering failure")
            logging.info(row)
    else:
        logging.info("  Skipped — provide --gff for per-locus assignment.")

    # ------------------------------------------------------------------ #
    # Check 3: Coverage distribution (5mC, valid_coverage)
    # ------------------------------------------------------------------ #
    section("CHECK 3: 5mC coverage distribution (valid_coverage)")
    for s, df in mc.items():
        cov = df["valid_coverage"]
        med, mean = cov.median(), cov.mean()
        p10, p90 = cov.quantile(0.10), cov.quantile(0.90)
        logging.info(f"  {s}: median={med:.1f}  mean={mean:.1f}  p10={p10:.1f}  p90={p90:.1f}")
        if med < 5:
            flag(f"{s}: median coverage {med:.1f} < 5")
        tsv_rows += [
            {"sample": s, "check": "cov_median", "value": round(med,  2)},
            {"sample": s, "check": "cov_mean",   "value": round(mean, 2)},
            {"sample": s, "check": "cov_p10",    "value": round(p10,  2)},
            {"sample": s, "check": "cov_p90",    "value": round(p90,  2)},
        ]

    # ------------------------------------------------------------------ #
    # Check 4: Methylation fraction distribution (5mC)
    # ------------------------------------------------------------------ #
    section("CHECK 4: 5mC methylation fraction (frac_modified)")
    for s, df in mc.items():
        mean_frac = df["frac_modified"].mean()
        logging.info(f"  {s}: mean frac_modified = {mean_frac:.4f}")
        if not (0.10 <= mean_frac <= 0.30):
            flag(f"{s}: mean frac_modified {mean_frac:.4f} outside expected range 0.10–0.30")
        tsv_rows.append({"sample": s, "check": "mean_frac_modified_5mC", "value": round(mean_frac, 4)})

    # ------------------------------------------------------------------ #
    # Check 5: Strand balance
    # ------------------------------------------------------------------ #
    section("CHECK 5: Strand balance (5mC sites)")
    for s, df in mc.items():
        counts = df["strand"].value_counts()
        n_plus  = int(counts.get("+", 0))
        n_minus = int(counts.get("-", 0))
        total = n_plus + n_minus
        pct_plus = 100 * n_plus / total if total else 0.0
        logging.info(f"  {s}: + = {n_plus:,} ({pct_plus:.1f}%)   - = {n_minus:,} ({100-pct_plus:.1f}%)")
        if total > 0 and not (40 <= pct_plus <= 60):
            flag(f"{s}: strand imbalance — {pct_plus:.1f}% plus-strand")
        tsv_rows.append({"sample": s, "check": "strand_pct_plus", "value": round(pct_plus, 1)})

    # ------------------------------------------------------------------ #
    # Check 6: Code composition (raw, before artefact exclusion)
    # ------------------------------------------------------------------ #
    section("CHECK 6: Modification code composition (raw, pre-exclusion)")
    for s, df in raw.items():
        counts = df["code"].value_counts()
        n_m   = int(counts.get(PRIMARY_CODE,   0))
        n_a   = int(counts.get(SECONDARY_CODE, 0))
        n_art = int(counts.get(ARTEFACT_CODE,  0))
        total = len(df)
        pct_a   = 100 * n_a   / total if total else 0.0
        pct_art = 100 * n_art / total if total else 0.0
        logging.info(f"  {s}: m={n_m:,}  a={n_a:,} ({pct_a:.2f}%)  21839={n_art:,} ({pct_art:.2f}%)")
        tsv_rows += [
            {"sample": s, "check": "raw_n_21839",   "value": n_art},
            {"sample": s, "check": "raw_pct_21839", "value": round(pct_art, 3)},
            {"sample": s, "check": "raw_pct_6mA",   "value": round(pct_a,   3)},
        ]

    # ------------------------------------------------------------------ #
    # Check 7: Cross-sample 5mC position overlap
    # ------------------------------------------------------------------ #
    section("CHECK 7: Cross-sample 5mC position overlap")
    pos_sets = {s: set(zip(df["chrom"], df["start"])) for s, df in mc.items()}
    sample_list = list(pos_sets)
    if len(sample_list) >= 2:
        shared = set.intersection(*pos_sets.values())
        union  = set.union(*pos_sets.values())
        pct_shared = 100 * len(shared) / len(union) if union else 0.0
        logging.info(f"  Shared across all {len(sample_list)} samples : {len(shared):,} sites")
        logging.info(f"  Union across all samples                  : {len(union):,} sites")
        logging.info(f"  Jaccard (shared / union)                  : {pct_shared:.1f}%")
        others = {s: set.union(*(v for k, v in pos_sets.items() if k != s)) for s in sample_list}
        for s in sample_list:
            private = len(pos_sets[s] - others[s])
            logging.info(f"  {s}: {len(pos_sets[s]):,} total,  {private:,} private")
            tsv_rows.append({"sample": s, "check": "private_5mC_sites", "value": private})
        if pct_shared < 20:
            flag(f"Very low cross-sample overlap ({pct_shared:.1f}%) — verify filtering used the same reference coordinates")
        tsv_rows += [
            {"sample": "all", "check": "shared_5mC_sites", "value": len(shared)},
            {"sample": "all", "check": "union_5mC_sites",  "value": len(union)},
            {"sample": "all", "check": "pct_shared",       "value": round(pct_shared, 1)},
        ]
    else:
        logging.info("  Skipped — need at least 2 samples for overlap analysis.")

    # ------------------------------------------------------------------ #
    # Summary + TSV output
    # ------------------------------------------------------------------ #
    section("SUMMARY")
    if WARNINGS:
        logging.info(f"  {len(WARNINGS)} warning(s):")
        for w in WARNINGS:
            logging.info(f"    - {w}")
    else:
        logging.info("  All checks passed — no warnings.")

    tsv_path = outdir / "sorgoleone_bedmethyl_validation.tsv"
    pd.DataFrame(tsv_rows).to_csv(tsv_path, sep="\t", index=False)
    logging.info(f"\n  Summary saved to: {tsv_path}")
    logging.info(f"  Log saved to:     {log_path}")


if __name__ == "__main__":
    main()
