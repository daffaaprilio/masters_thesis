#!/usr/bin/env python3
"""
Visualize methylation sequence context (CpG / CHG / CHH) from the output of
methylation_context_samples.sh.

Outputs (analysis/02_methylation_landscape/figures/):
  fig_context_site_proportion.png   – stacked bar: fraction of 5mC sites per context
  fig_context_methylation_level.png – grouped bar: mean % methylation per context

Usage:
    ./docker/run.sh python3 workflow/scripts/methylation_context_viz.py [path/to/tsv]
    # Omit path to use the most recent TSV in analysis/data/methylation_context/
"""

import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from datetime import datetime
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR  = REPO_ROOT / "analysis/data/methylation_context"
OUT_DIR   = REPO_ROOT / "analysis/02_methylation_landscape/figures"
LOG_DIR   = REPO_ROOT / "analysis/logs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
_log_path  = LOG_DIR / f"methylation_context_viz_{_timestamp}.log"
_root = logging.getLogger()
_root.setLevel(logging.INFO)
_ch = logging.StreamHandler()
_ch.setFormatter(logging.Formatter("%(message)s"))
_root.addHandler(_ch)
_fh = logging.FileHandler(_log_path)
_fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
_root.addHandler(_fh)
logging.info(f"Log: {_log_path}")

SAMPLES  = ["SBC4", "SBC10", "SBC11", "SBC23"]
CONTEXTS = ["CpG", "CHG", "CHH"]

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":          9,
    "axes.labelsize":     9,
    "axes.titlesize":    10,
    "legend.fontsize":    7.5,
    "figure.dpi":       150,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

CTX_COLORS = {"CpG": "#2166ac", "CHG": "#4dac26", "CHH": "#d01c8b"}

# ── Load ───────────────────────────────────────────────────────────────────────

def find_tsv() -> Path:
    tsvs = sorted(DATA_DIR.glob("methylation_context_*.tsv"))
    if not tsvs:
        raise FileNotFoundError(f"No methylation_context TSV found in {DATA_DIR}")
    return tsvs[-1]


def save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved {name}")


tsv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else find_tsv()
logging.info(f"Reading: {tsv_path}")

df = pd.read_csv(tsv_path, sep="\t")
df["sample"]  = pd.Categorical(df["sample"],  categories=SAMPLES,  ordered=True)
df["context"] = pd.Categorical(df["context"], categories=CONTEXTS, ordered=True)
df = df.sort_values(["sample", "context"])

present = [s for s in SAMPLES if s in df["sample"].values]

# ── fig_context_site_proportion – site proportion per context (stacked bar) ──

fig, ax = plt.subplots(figsize=(6, 4.5))

pivot = (df.pivot_table(index="sample", columns="context",
                        values="n_sites", fill_value=0)
           .reindex(index=present)
           .reindex(columns=CONTEXTS, fill_value=0))

prop = pivot.div(pivot.sum(axis=1), axis=0)

bottom = np.zeros(len(prop))
for ctx in CONTEXTS:
    vals = prop[ctx].values
    ax.bar(range(len(prop)), vals, bottom=bottom,
           color=CTX_COLORS[ctx], label=ctx, width=0.6,
           edgecolor="white", linewidth=0.5)
    for i, (v, b) in enumerate(zip(vals, bottom)):
        if v > 0.02:
            ax.text(i, b + v / 2, f"{v:.1%}",
                    ha="center", va="center",
                    fontsize=7.5, color="white", fontweight="bold")
    bottom += vals

ax.set_xticks(range(len(prop)))
ax.set_xticklabels(prop.index, fontsize=8)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax.set_ylabel("Proportion of 5mC sites")
ax.set_title("5mC site distribution by sequence context", fontweight="bold")
ax.legend(title="Context", fontsize=8, title_fontsize=8,
          loc="upper right", framealpha=0.6)

for i, s in enumerate(prop.index):
    total = int(pivot.loc[s].sum())
    ax.text(i, -0.06, f"n={total:,}", ha="center", va="top",
            fontsize=6.5, color="#555555", transform=ax.get_xaxis_transform())

plt.tight_layout()
save(fig, "fig_context_site_proportion.png")

# ── fig_context_methylation_level – mean methylation level per context ───────

fig, ax = plt.subplots(figsize=(7, 4.5))

SAMPLE_PAL = dict(zip(SAMPLES, sns.color_palette("Set2", len(SAMPLES))))

n_ctx     = len(CONTEXTS)
n_samples = len(present)
group_w   = 0.8
bar_w     = group_w / n_samples
offsets   = np.linspace(-(group_w - bar_w) / 2, (group_w - bar_w) / 2, n_samples)

for j, sample in enumerate(present):
    sdf = df[df["sample"] == sample].set_index("context")
    vals = [sdf.loc[ctx, "mean_pct_mod"] if ctx in sdf.index else 0
            for ctx in CONTEXTS]
    xs = np.arange(n_ctx) + offsets[j]
    ax.bar(xs, vals, width=bar_w, color=SAMPLE_PAL[sample],
           label=sample, edgecolor="white", linewidth=0.4)

ax.set_xticks(np.arange(n_ctx))
ax.set_xticklabels(CONTEXTS, fontsize=9)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=1))
ax.set_ylabel("Mean methylation level (%)")
ax.set_title("Mean 5mC methylation level per sequence context", fontweight="bold")
ax.legend(fontsize=7.5, framealpha=0.6, ncol=2)

plt.tight_layout()
save(fig, "fig_context_methylation_level.png")

logging.info(f"All figures saved -> {OUT_DIR}/")
