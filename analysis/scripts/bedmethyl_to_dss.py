#!/usr/bin/env python3
"""
Convert filtered sorgoleone bedMethyl files to DSS input format.

For 5mC: collapses Watson/Crick CpG strand pairs into single positions
by summing N (valid_coverage) and X (n_mod) from both strands.
For 6mA: keeps strands separate (6mA is not symmetric).

Output: {sample}.5mC.dss.txt and {sample}.6mA.dss.txt

Usage:
    python analysis/scripts/bedmethyl_to_dss.py \
        [--indir  analysis/data/sorgoleone_bedmethyl] \
        [--outdir analysis/data/sorgoleone_DSS]
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

ARTEFACT_CODE = "21839"
PRIMARY_CODE  = "m"   # 5mC
SECONDARY_CODE = "a"  # 6mA
MIN_COVERAGE  = 5

DSS_COLS = ["chr", "pos", "N", "X"]


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


def load_sample(path):
    return pd.read_csv(
        path, sep="\t", header=None, names=BEDMETHYL_COLS,
        dtype={"code": str, "valid_coverage": float, "n_mod": float,
               "frac_modified": float},
    )


def check_intervals(df, sample):
    multi = df[df["end"] - df["start"] != 1]
    if not multi.empty:
        logging.warning(f"*** {len(multi):,} multi-base intervals in {sample} — review expected")


def check_counts(dss, label):
    bad = dss[dss["X"] > dss["N"]]
    if not bad.empty:
        logging.warning(f"*** {len(bad):,} rows where X > N in {label}")
    neg = dss[(dss["N"] < 0) | (dss["X"] < 0)]
    if not neg.empty:
        logging.warning(f"*** {len(neg):,} rows with negative N or X in {label}")


def check_duplicates(dss, label):
    dupes = dss.duplicated(subset=["chr", "pos"])
    if dupes.any():
        logging.warning(f"*** {dupes.sum():,} duplicate positions in {label}")


def collapse_cpg_strands(df):
    """
    Collapse Watson (+) and Crick (-) strand rows for the same CpG.

    In modkit pileup output (BED, 0-based):
      + strand C at position X   → start = X
      - strand C at position X+1 → start = X+1  (complement of G at X on + strand)

    Canonical DSS position (1-based) for both = X+1:
      + strand: pos = start + 1
      - strand: pos = start      (start - 1 is 0-based partner C; +1 gives start)

    Sum N and X across both strands at each collapsed position.
    """
    plus  = df[df["strand"] == "+"].copy()
    minus = df[df["strand"] == "-"].copy()

    plus ["pos"] = plus ["start"] + 1   # 0-based → 1-based
    minus["pos"] = minus["start"]        # already aligns to + strand 1-based position

    combined = pd.concat([
        plus [["chrom", "pos", "valid_coverage", "n_mod"]],
        minus[["chrom", "pos", "valid_coverage", "n_mod"]],
    ])

    collapsed = (
        combined
        .groupby(["chrom", "pos"], sort=False)[["valid_coverage", "n_mod"]]
        .sum()
        .reset_index()
        .rename(columns={"chrom": "chr", "valid_coverage": "N", "n_mod": "X"})
    )
    collapsed["N"] = collapsed["N"].astype(int)
    collapsed["X"] = collapsed["X"].astype(int)
    return collapsed[DSS_COLS].sort_values(["chr", "pos"]).reset_index(drop=True)


def to_dss_no_collapse(df):
    """Convert to DSS format without strand collapsing (for 6mA)."""
    out = pd.DataFrame({
        "chr": df["chrom"],
        "pos": (df["start"] + 1).astype(int),
        "N":   df["valid_coverage"].astype(int),
        "X":   df["n_mod"].astype(int),
    })
    return out.sort_values(["chr", "pos"]).reset_index(drop=True)


def section(title):
    logging.info(f"\n{'─' * 60}\n  {title}\n{'─' * 60}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indir",  default="analysis/data/sorgoleone_bedmethyl",
                        help="Directory with .sorgoleone.bed files")
    parser.add_argument("--outdir", default="analysis/data/sorgoleone_DSS",
                        help="Output directory for DSS files")
    parser.add_argument("--min-coverage", type=int, default=MIN_COVERAGE,
                        help=f"Minimum valid_coverage to keep (default: {MIN_COVERAGE})")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = outdir / f"bedmethyl_to_dss_{timestamp}.log"
    setup_logging(log_path)

    logging.info(f"Log: {log_path}")
    logging.info(f"Input dir: {args.indir}")
    logging.info(f"Min coverage filter: {args.min_coverage}")

    indir = Path(args.indir)
    bed_files = sorted(indir.glob("*.sorgoleone.bed"))
    if not bed_files:
        logging.error(f"No .sorgoleone.bed files found in {indir}")
        sys.exit(1)

    samples = {p.name.split(".")[0]: p for p in bed_files}
    logging.info(f"Samples: {', '.join(samples)}")

    summary_rows = []

    for sample, path in samples.items():
        section(sample)

        raw = load_sample(path)
        check_intervals(raw, sample)

        clean = raw[raw["code"] != ARTEFACT_CODE]

        for code, label, do_collapse in [
            (PRIMARY_CODE,   "5mC", True),
            (SECONDARY_CODE, "6mA", False),
        ]:
            ctx = clean[clean["code"] == code].copy()
            n_raw = len(ctx)

            if n_raw == 0:
                logging.info(f"  {label}: no sites — skipping")
                continue

            ctx = ctx[ctx["valid_coverage"] >= args.min_coverage]
            n_after_cov = len(ctx)
            logging.info(f"  {label}: {n_raw:,} raw → {n_after_cov:,} after coverage ≥ {args.min_coverage} filter")

            if n_after_cov == 0:
                logging.info(f"  {label}: no sites remain after coverage filter — skipping")
                continue

            if do_collapse:
                n_before_collapse = n_after_cov
                dss = collapse_cpg_strands(ctx)
                n_after_collapse = len(dss)
                logging.info(f"  {label}: {n_before_collapse:,} strand rows → {n_after_collapse:,} after CpG collapse")
            else:
                dss = to_dss_no_collapse(ctx)
                n_after_collapse = len(dss)

            check_counts(dss, f"{sample} {label}")
            check_duplicates(dss, f"{sample} {label}")

            out_path = outdir / f"{sample}.{label}.dss.txt"
            dss.to_csv(out_path, sep="\t", index=False)
            logging.info(f"  {label}: wrote {out_path.name}")

            mean_frac = (dss["X"] / dss["N"]).mean()
            summary_rows.append({
                "sample":                      sample,
                "context":                     label,
                "sites_raw":                   n_raw,
                "sites_after_coverage_filter": n_after_cov,
                "sites_after_collapse":        n_after_collapse,
                "min_N":                       int(dss["N"].min()),
                "median_N":                    dss["N"].median(),
                "max_N":                       int(dss["N"].max()),
                "mean_frac":                   round(mean_frac, 4),
            })

    summary_path = outdir / "conversion_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    logging.info(f"\nConversion summary saved to: {summary_path}")
    logging.info(f"Log saved to: {log_path}")


if __name__ == "__main__":
    main()
