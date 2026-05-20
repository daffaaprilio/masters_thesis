#!/usr/bin/env python3
"""
Visualize modkit summary TSV output.

Produces one figure in analysis/02_methylation_landscape/figures/:
  fig_modkit_summary_composition.png – stacked modification composition per base

Usage:
    ./docker/run.sh python3 analysis/scripts/modkit_summary_viz.py [path/to/tsv]
    # Omit path to use the most recent TSV in analysis/data/modkit_summary/
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR  = REPO_ROOT / "analysis/data/modkit_summary"
OUT_DIR   = REPO_ROOT / "analysis/02_methylation_landscape/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Sample metadata ────────────────────────────────────────────────────────────
SAMPLES   = ["SBC4", "SBC10", "SBC11", "SBC23"]

# ── Style ──────────────────────────────────────────────────────────────────────
CODE_COLORS = {
    "m":     "#e07b54",
    "21839": "#b59fd4",
    "-":     "#d0d0d0",
    "a":     "#5bafd4",
}
CODE_LABELS = {"m": "5mC", "a": "6mA", "-": "unmod", "21839": "other"}

plt.rcParams.update({
    "figure.dpi":          150,
    "font.size":           9,
    "axes.spines.top":     False,
    "axes.spines.right":   False,
})


def find_tsv() -> Path:
    tsvs = sorted(DATA_DIR.glob("modkit_summary_*.tsv"))
    if not tsvs:
        raise FileNotFoundError(f"No modkit_summary TSV found in {DATA_DIR}")
    return tsvs[-1]


def save(fig: plt.Figure, name: str) -> None:
    out = OUT_DIR / name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close(fig)


# ── Load ───────────────────────────────────────────────────────────────────────
tsv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else find_tsv()
print(f"Reading: {tsv_path}")

df = pd.read_csv(tsv_path, sep="\t", dtype={"code": str})
df["sample"] = pd.Categorical(df["sample"], categories=SAMPLES, ordered=True)
df = df.sort_values("sample")

present = [s for s in SAMPLES if s in df["sample"].values]

# ── Modification composition (stacked) ────────────────────────────────────────
fig, (ax_c, ax_a) = plt.subplots(1, 2, figsize=(10, 4))

for ax, base, title in [
    (ax_c, "C", "Cytosine call composition"),
    (ax_a, "A", "Adenine call composition"),
]:
    sub = df[df["base"] == base].copy()
    sub = sub[sub["sample"].isin(present)]

    # Code order: modifications first, unmodified last
    codes = [c for c in ["m", "a", "21839", "-"] if c in sub["code"].values]

    pivot = (sub.pivot_table(index="sample", columns="code",
                             values="pass_frac", fill_value=0)
               .reindex(index=present)
               .reindex(columns=codes, fill_value=0))

    bottom = np.zeros(len(pivot))
    for code in codes:
        ax.bar(range(len(pivot)), pivot[code].values, bottom=bottom,
               color=CODE_COLORS.get(code, "#aaaaaa"),
               label=CODE_LABELS.get(code, code),
               width=0.6, edgecolor="white", linewidth=0.5)
        # Label segments > 1 %
        for i, (v, b) in enumerate(zip(pivot[code].values, bottom)):
            if v > 0.01:
                ax.text(i, b + v / 2, f"{v:.1%}",
                        ha="center", va="center", fontsize=7,
                        color="white", fontweight="bold")
        bottom += pivot[code].values

    ax.set_xticks(range(len(pivot)))
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_ylabel("Fraction of pass-threshold calls")
    ax.set_title(title, fontweight="bold")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.5)

fig.suptitle("Modification composition per base (pass-threshold calls)", y=1.01, fontsize=11)
fig.tight_layout()
save(fig, "fig_modkit_summary_composition.png")
