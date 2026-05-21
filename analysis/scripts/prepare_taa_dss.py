#!/usr/bin/env python3
"""
Convert whole-genome filtered bedMethyl files to DSS input format for TAA DMR analysis.

Reads resources/bedmethyl/{sample}.filtered.bed (already filtered for ≥10
valid reads by modkit) for each of the four sorghum accessions, collapses
Watson/Crick CpG strand pairs for 5mC, and writes DSS-ready text files.

Run via:
    ./docker/run.sh python3 analysis/scripts/prepare_taa_dss.py

Output: analysis/data/taa_DSS/{sample}.5mC.dss.txt
        analysis/data/taa_DSS/{sample}.6mA.dss.txt
        analysis/data/taa_DSS/conversion_summary.tsv
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# ── Config ─────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).resolve().parent.parent.parent
BED_DIR   = ROOT / "resources/bedmethyl"
OUT_DIR   = ROOT / "analysis/data/taa_DSS"
LOG_DIR   = ROOT / "analysis/logs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]

BEDMETHYL_COLS = [
    "chrom", "start", "end", "code", "coverage", "strand",
    "start2", "end2", "color", "valid_coverage", "frac_modified",
    "n_mod", "n_canonical", "n_other_mod", "n_delete",
    "n_fail", "n_diff", "n_nocall",
]
ARTEFACT_CODE  = "21839"
CODE_5MC       = "m"
CODE_6MA       = "a"
DSS_COLS       = ["chr", "pos", "N", "X"]

_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(LOG_DIR / f"prepare_taa_dss_{_ts}.log"),
    ],
)


def load_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path, sep="\t", header=None, names=BEDMETHYL_COLS,
        dtype={"code": str, "valid_coverage": float, "n_mod": float,
               "frac_modified": float},
    )


def collapse_cpg_strands(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse Watson (+) and Crick (−) CpG pairs into a single 1-based position.

    Modkit BED convention (0-based):
      + strand C at position X → start = X  → DSS pos = X + 1
      − strand C at position X+1 → start = X+1 → DSS pos = X+1 (aligns to + strand)

    N and X are summed across both strands.
    """
    plus  = df[df["strand"] == "+"].copy()
    minus = df[df["strand"] == "-"].copy()
    plus ["pos"] = plus ["start"] + 1
    minus["pos"] = minus["start"]

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


def to_dss_simple(df: pd.DataFrame) -> pd.DataFrame:
    """DSS conversion without strand collapsing (for 6mA)."""
    out = pd.DataFrame({
        "chr": df["chrom"],
        "pos": (df["start"] + 1).astype(int),
        "N":   df["valid_coverage"].astype(int),
        "X":   df["n_mod"].astype(int),
    })
    return out.sort_values(["chr", "pos"]).reset_index(drop=True)


def check(dss: pd.DataFrame, label: str) -> None:
    bad = dss[dss["X"] > dss["N"]]
    if not bad.empty:
        logging.warning(f"  *** {len(bad):,} rows where X > N in {label}")
    dupes = dss.duplicated(subset=["chr", "pos"])
    if dupes.any():
        logging.warning(f"  *** {dupes.sum():,} duplicate positions in {label}")


# ── Main ───────────────────────────────────────────────────────────────────────
logging.info(f"Input:  {BED_DIR}")
logging.info(f"Output: {OUT_DIR}")

summary_rows = []

for sample in SAMPLES:
    bed_path = BED_DIR / f"{sample}.filtered.bed"
    if not bed_path.exists():
        logging.error(f"Missing: {bed_path}")
        sys.exit(1)

    logging.info(f"\n{'─' * 60}\n  {sample}\n{'─' * 60}")
    raw = load_bed(bed_path)
    logging.info(f"  Loaded {len(raw):,} rows")

    clean = raw[raw["code"] != ARTEFACT_CODE]

    for code, label, collapse in [
        (CODE_5MC, "5mC", True),
        (CODE_6MA, "6mA", False),
    ]:
        ctx = clean[clean["code"] == code].copy()
        if ctx.empty:
            logging.info(f"  {label}: no rows — skipping")
            continue

        logging.info(f"  {label}: {len(ctx):,} rows")

        dss = collapse_cpg_strands(ctx) if collapse else to_dss_simple(ctx)
        if collapse:
            logging.info(f"  {label}: {len(ctx):,} strand rows → {len(dss):,} CpG positions after collapse")

        check(dss, f"{sample} {label}")

        out_path = OUT_DIR / f"{sample}.{label}.dss.txt"
        dss.to_csv(out_path, sep="\t", index=False)
        logging.info(f"  {label}: wrote {out_path.name}  ({len(dss):,} sites)")

        summary_rows.append({
            "sample":          sample,
            "context":         label,
            "rows_in":         len(ctx),
            "sites_out":       len(dss),
            "mean_coverage":   round(dss["N"].mean(), 1),
            "mean_frac_mod":   round((dss["X"] / dss["N"]).mean(), 4),
        })

summary_path = OUT_DIR / "conversion_summary.tsv"
pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
logging.info(f"\nConversion summary → {summary_path}")
logging.info("Done.")
