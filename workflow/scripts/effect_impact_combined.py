#!/usr/bin/env python3
"""
Effect Impact Distribution, computed directly from a combined multi-sample
SnpEff-annotated VCF (results/combined/all.annotated.vcf.gz), rather than
from per-sample SnpEff csvStats files (see variant_landscape.py fig02).

For each variant, the worst (highest-severity) impact across all ANN entries
is taken as that variant's impact. A sample is counted as carrying a variant
if its genotype contains at least one ALT allele (not missing, not 0/0).

Output: analysis/figures/fig02_effect_impact_combined.png
"""

import argparse
import logging
from pathlib import Path

import cyvcf2
import numpy as np
import pandas as pd

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]

IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
IMPACT_RANK = {imp: i for i, imp in enumerate(IMPACT_ORDER)}  # lower = worse
IMPACT_COLORS = {
    "HIGH":     "#c44e52",
    "MODERATE": "#dd8452",
    "LOW":      "#4c72b0",
    "MODIFIER": "#8172b2",
}


def worst_impact(ann_raw: str) -> str | None:
    worst = None
    for annotation in ann_raw.split(","):
        parts = annotation.split("|")
        if len(parts) < 3:
            continue
        impact = parts[2]
        rank = IMPACT_RANK.get(impact)
        if rank is None:
            continue
        if worst is None or rank < IMPACT_RANK[worst]:
            worst = impact
    return worst


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--vcf",
        default="results/combined/all.annotated.vcf.gz",
        help="Combined SnpEff-annotated multi-sample VCF",
    )
    parser.add_argument(
        "-o", "--out",
        default="analysis/figures/fig02_effect_impact_combined.png",
        help="Output figure path",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger(__name__)

    vcf = cyvcf2.VCF(args.vcf)
    samples = vcf.samples
    log.info("Samples in VCF: %s", samples)
    sample_idx = {s: samples.index(s) for s in SAMPLES if s in samples}

    counts = {s: {imp: 0 for imp in IMPACT_ORDER} for s in sample_idx}

    n = 0
    for variant in vcf:
        ann_raw = variant.INFO.get("ANN")
        if not ann_raw:
            continue
        impact = worst_impact(ann_raw)
        if impact is None:
            continue

        genotypes = variant.genotypes
        for s, idx in sample_idx.items():
            a1, a2, _phased = genotypes[idx]
            if a1 > 0 or a2 > 0:
                counts[s][impact] += 1

        n += 1
        if n % 500_000 == 0:
            log.info("Processed %d variants", n)

    log.info("Total variants processed: %d", n)

    df_imp = pd.DataFrame(counts).T[IMPACT_ORDER]
    df_imp = df_imp.reindex(SAMPLES).dropna(how="all")
    pct_imp = df_imp.div(df_imp.sum(axis=1), axis=0) * 100

    log.info("Impact counts:\n%s", df_imp)
    log.info("Impact percentages:\n%s", pct_imp.round(3))

    plot_effect_impact(pct_imp, Path(args.out))
    log.info("Saved figure to %s", args.out)


def plot_effect_impact(df_imp: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family":     "DejaVu Sans",
        "font.size":       9,
        "axes.labelsize":  9,
        "axes.titlesize":  10,
        "legend.fontsize": 7.5,
        "figure.dpi":      150,
    })

    max_small = (df_imp.drop(columns="MODIFIER", errors="ignore")
                       .sum(axis=1).max())
    LO_MAX = np.ceil(max_small * 1.35)
    HI_MIN = 93.0

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, figsize=(6, 5),
        gridspec_kw={"height_ratios": [1, 3], "hspace": 0.07},
    )

    x = np.arange(len(df_imp))
    for ax in (ax_top, ax_bot):
        bottom = np.zeros(len(df_imp))
        for col in df_imp.columns:
            ax.bar(x, df_imp[col].values, bottom=bottom,
                   color=IMPACT_COLORS.get(col, "#aaaaaa"), label=col, width=0.6)
            bottom += df_imp[col].values

    ax_bot.set_ylim(0, LO_MAX)
    ax_top.set_ylim(HI_MIN, 101)

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(df_imp.index, rotation=30, ha="right")
    ax_top.tick_params(axis="x", bottom=False)

    ax_top.spines.bottom.set_visible(False)
    ax_bot.spines.top.set_visible(False)

    d = 0.018
    kw = dict(color="k", clip_on=False, lw=1.2)
    ax_top.plot([-d, +d], [-d, +d], transform=ax_top.transAxes, **kw)
    ax_top.plot([1 - d, 1 + d], [-d, +d], transform=ax_top.transAxes, **kw)
    ax_bot.plot([-d, +d], [1 - d, 1 + d], transform=ax_bot.transAxes, **kw)
    ax_bot.plot([1 - d, 1 + d], [1 - d, 1 + d], transform=ax_bot.transAxes, **kw)

    ax_bot.set_ylabel("% of effects")
    ax_top.set_title("Effect Impact Distribution")
    ax_top.legend(loc="lower right", framealpha=0.8, fontsize=7.5)

    plt.tight_layout()

    bot_top = ax_bot.get_position().y1
    top_bot = ax_top.get_position().y0
    x_left = ax_bot.get_position().x0
    fig.text(x_left - 0.04, (bot_top + top_bot) / 2, "⋮",
             ha="center", va="center", fontsize=13, color="gray",
             transform=fig.transFigure, clip_on=False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
