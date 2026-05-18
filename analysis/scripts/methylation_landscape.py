#!/usr/bin/env python3
"""
Methylation landscape EDA from modkit bedMethyl files.

Input:  resources/bedmethyl/{sample}.filtered.bed  (falls back to {sample}.bed)
Output: analysis/02_methylation_landscape/figures/
  fig01_methyl_distribution.png  – per-sample per-context pct_mod KDE
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from scipy.stats import gaussian_kde

# ── Config ──────────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).parent.parent.parent / "resources/bedmethyl"
OUT_DIR  = Path(__file__).parent.parent / "02_methylation_landscape/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES  = ["SBC4", "SBC10", "SBC11", "SBC23"]
MIN_COV  = 10   # minimum N_valid reads per cytosine position

BEDMETHYL_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thick_start", "thick_end", "color",
    "N_valid", "pct_mod", "N_mod", "N_canonical",
    "N_nocall", "N_filtered", "N_diff", "N_delete",
]

CONTEXT_LABELS = {"m": "5mC", "h": "5hmC", "a": "6mA"}

SAMPLE_COLORS = dict(zip(SAMPLES, sns.color_palette("Set2", len(SAMPLES))))

plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       9,
    "axes.labelsize":  9,
    "axes.titlesize":  10,
    "legend.fontsize": 7.5,
    "figure.dpi":      150,
})

# ── Data loading ─────────────────────────────────────────────────────────────

def load_bedmethyl(sample: str) -> pd.DataFrame:
    for fname in (f"{sample}.filtered.bed", f"{sample}.bed"):
        path = DATA_DIR / fname
        if path.exists():
            df = pd.read_csv(
                path, sep="\t", header=None, names=BEDMETHYL_COLS,
                usecols=["chrom", "name", "N_valid", "pct_mod"],
                dtype={"chrom": str, "name": str, "N_valid": int, "pct_mod": float},
            )
            df = df[df["N_valid"] >= MIN_COV].copy()
            df["sample"] = sample
            return df
    print(f"  [warn] {sample}: no bedMethyl file in {DATA_DIR}", file=sys.stderr)
    return pd.DataFrame()


data = {s: load_bedmethyl(s) for s in SAMPLES}
available = {s: df for s, df in data.items() if not df.empty}

if not available:
    print(
        f"No bedMethyl data found under {DATA_DIR}.\n"
        "Run the methylation calling pipeline first:\n"
        "  ./docker/snakemake.sh methylation_all --cores 24"
    )
    sys.exit(0)

combined  = pd.concat(available.values(), ignore_index=True)
contexts  = sorted(combined["name"].unique())
n_ctx     = len(contexts)

# ── Helpers ───────────────────────────────────────────────────────────────────

def save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(name)


def kde_curve(values: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    if len(values) < 2:
        return np.zeros_like(x_grid)
    bw = max(0.5, np.std(values) * len(values) ** -0.2)
    return gaussian_kde(values, bw_method=bw / np.std(values))(x_grid)

# ── fig01 – per-sample per-context methylation distribution ──────────────────
# One column per modification context; KDE of pct_mod per sample overlaid.

x_grid = np.linspace(0, 100, 500)

fig, axes = plt.subplots(
    1, n_ctx,
    figsize=(4.5 * n_ctx, 4.2),
    sharey=False,
    squeeze=False,
)

for col, ctx in enumerate(contexts):
    ax = axes[0, col]
    label = CONTEXT_LABELS.get(ctx, ctx)

    for sample, df in available.items():
        vals = df.loc[df["name"] == ctx, "pct_mod"].values
        if len(vals) == 0:
            continue
        y = kde_curve(vals, x_grid)
        color = SAMPLE_COLORS[sample]
        ax.plot(x_grid, y, lw=1.8, color=color, label=sample)
        ax.fill_between(x_grid, y, alpha=0.12, color=color)

    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Methylation (%)")
    ax.set_title(f"{label} context")
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
    ax.tick_params(axis="x", which="minor", length=3)

    n_sites = {
        s: (available[s]["name"] == ctx).sum()
        for s in available
    }
    note = "\n".join(f"{s}: {n:,}" for s, n in n_sites.items() if n > 0)
    ax.text(
        0.97, 0.97, note,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=6.5, color="#555555",
        family="monospace",
    )

    if col == 0:
        ax.set_ylabel("Density")

# shared legend on the right
legend_handles = [
    mpatches.Patch(color=SAMPLE_COLORS[s], label=s)
    for s in available
]
axes[0, -1].legend(
    handles=legend_handles,
    loc="upper left",
    framealpha=0.8,
    title="Sample",
)

cov_note = f"Coverage filter: N_valid ≥ {MIN_COV}"
fig.text(0.5, -0.02, cov_note, ha="center", fontsize=7.5, color="#888888")
fig.suptitle("Per-Sample Methylation Level Distribution", y=1.02, fontsize=11)

plt.tight_layout()
save(fig, "fig01_methyl_distribution.png")

print(f"\nAll figures saved → {OUT_DIR}/")
