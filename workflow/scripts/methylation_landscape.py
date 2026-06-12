#!/usr/bin/env python3
"""
Methylation landscape EDA from modkit bedMethyl files.

Input:  resources/bedmethyl/{sample}.filtered.bed  (falls back to {sample}.bed)
Output: analysis/02_methylation_landscape/figures/
  fig01_context_composition.png      – stacked bar: proportion of modification contexts per sample
  fig02_5mc_by_seq_context.png       – violin: 5mC frequency distribution per CpG/CHG/CHH per sample
  fig03_chrplot_{CpG,CHG,CHH}.png    – chromosome-level methylation density per sequence context
"""

import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from pathlib import Path
from pyfaidx import Fasta

# ── Config ──────────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).parent.parent.parent / "resources/bedmethyl"
OUT_DIR  = Path(__file__).parent.parent / "02_methylation_landscape/figures"
LOG_DIR  = Path(__file__).parent.parent / "logs"
REF_FA   = Path(__file__).parent.parent.parent / "resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
_log_path  = LOG_DIR / f"methylation_landscape_{_timestamp}.log"
_root = logging.getLogger()
_root.setLevel(logging.INFO)
logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().handlers[0].setFormatter(logging.Formatter("%(message)s"))
_fh = logging.FileHandler(_log_path)
_fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
_root.addHandler(_fh)
logging.info(f"Log: {_log_path}")

SAMPLES  = ["SBC4", "SBC10", "SBC11", "SBC23"]
MIN_COV  = 10   # minimum N_valid reads per cytosine position

BEDMETHYL_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thick_start", "thick_end", "color",
    "N_valid", "pct_mod", "N_mod", "N_canonical",
    "N_other_mod", "N_delete", "N_filtered", "N_diff", "N_nocall",
]

CONTEXT_LABELS = {"m": "5mC", "h": "5hmC", "a": "6mA"}

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
                usecols=["chrom", "start", "strand", "name", "N_valid", "N_mod", "pct_mod"],
                dtype={"chrom": str, "start": int, "strand": str, "name": str, "N_valid": int, "N_mod": int, "pct_mod": float},
            )
            df = df[df["N_valid"] >= MIN_COV].copy()
            df["sample"] = sample
            return df
    logging.warning(f"{sample}: no bedMethyl file in {DATA_DIR}")
    return pd.DataFrame()


data = {s: load_bedmethyl(s) for s in SAMPLES}
available = {s: df for s, df in data.items() if not df.empty}

if not available:
    logging.error(
        f"No bedMethyl data found under {DATA_DIR}. "
        "Run the methylation calling pipeline first: "
        "./docker/snakemake.sh methylation_all --cores 24"
    )
    sys.exit(0)

combined = pd.concat(available.values(), ignore_index=True)
contexts = sorted(combined["name"].unique())

# ── Helpers ───────────────────────────────────────────────────────────────────

# Complement lookup table (ASCII index → ASCII complement)
_COMP = np.zeros(256, dtype=np.uint8)
for _a, _b in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    _COMP[_a] = _b


def assign_seq_context(df: pd.DataFrame, fa: Fasta) -> pd.Series:
    """Classify each 5mC position as CpG, CHG, or CHH using the reference."""
    G = ord("G")
    result = pd.Series(index=df.index, dtype=str)
    for chrom, grp in df.groupby("chrom"):
        try:
            seq_b = np.frombuffer(fa[chrom][:].upper().encode(), dtype=np.uint8)
        except KeyError:
            continue
        n   = len(seq_b)
        pos = grp["start"].values.astype(int)
        is_plus = grp["strand"].values == "+"

        b1 = np.full(len(grp), ord("N"), dtype=np.uint8)
        b2 = np.full(len(grp), ord("N"), dtype=np.uint8)

        pi = np.where(is_plus)[0]
        if pi.size:
            p = pos[pi]
            m1, m2 = p + 1 < n, p + 2 < n
            b1[pi[m1]] = seq_b[p[m1] + 1]
            b2[pi[m2]] = seq_b[p[m2] + 2]

        mi = np.where(~is_plus)[0]
        if mi.size:
            p = pos[mi]
            m1, m2 = p - 1 >= 0, p - 2 >= 0
            b1[mi[m1]] = _COMP[seq_b[p[m1] - 1]]
            b2[mi[m2]] = _COMP[seq_b[p[m2] - 2]]

        result.loc[grp.index] = np.where(b1 == G, "CpG", np.where(b2 == G, "CHG", "CHH"))
    return result


def save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved {name}")


# ── fig01 – methylation context composition per sample ───────────────────────
# Stacked bar: proportion of total N_mod attributed to each modification context.

CTX_COLORS = dict(zip(contexts, sns.color_palette("Set1", len(contexts))))

ctx_totals = (
    combined.groupby(["sample", "name"])["N_mod"]
    .sum()
    .unstack(fill_value=0)
)
ctx_props = ctx_totals.div(ctx_totals.sum(axis=1), axis=0)
ctx_props = ctx_props.reindex([s for s in SAMPLES if s in available])

fig, ax = plt.subplots(figsize=(max(4.5, 1.4 * len(available)), 4.2))

bottom = np.zeros(len(ctx_props))
for ctx in ctx_props.columns:
    label = CONTEXT_LABELS.get(ctx, ctx)
    vals  = ctx_props[ctx].values
    bars  = ax.bar(ctx_props.index, vals, bottom=bottom,
                   color=CTX_COLORS[ctx], label=label, width=0.55)
    for bar, v in zip(bars, vals):
        if v >= 0.03:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_y() + v / 2,
                f"{v:.1%}",
                ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold",
            )
    bottom += vals

ax.set_ylim(0, 1)
ax.set_ylabel("Proportion of modified base calls")
ax.set_xlabel("Sample")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:.0%}"))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.legend(title="Context", bbox_to_anchor=(1.01, 1), loc="upper left", framealpha=0.8)
ax.set_title("Methylation Context Composition per Sample")

cov_note = f"Coverage filter: N_valid ≥ {MIN_COV}"
fig.text(0.5, -0.02, cov_note, ha="center", fontsize=7.5, color="#888888")

plt.tight_layout()
save(fig, "fig01_context_composition.png")

# ── fig02 – 5mC frequency distribution by sequence context (CpG/CHG/CHH) ─────
# Violin plot of methylation frequency (0–1) per context per sample.

if not REF_FA.exists():
    logging.warning(f"reference not found at {REF_FA}; skipping fig02")
else:
    mc = combined[combined["name"] == "m"].copy()
    mc["freq"] = mc["pct_mod"] / 100.0

    logging.info("assigning sequence contexts from reference …")
    fa = Fasta(str(REF_FA), as_raw=True, sequence_always_upper=True)
    mc["seq_ctx"] = assign_seq_context(mc, fa)

    SEQ_CTXS     = ["CpG", "CHG", "CHH"]
    SAMPLE_ORDER = [s for s in SAMPLES if s in available]
    SAMPLE_PAL   = dict(zip(SAMPLES, sns.color_palette("Set2", len(SAMPLES))))
    MAX_PTS      = 60_000  # cap per context for KDE speed

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True)

    for ax, ctx in zip(axes, SEQ_CTXS):
        sub = mc[mc["seq_ctx"] == ctx]
        if len(sub) > MAX_PTS:
            sub = sub.sample(MAX_PTS, random_state=42)

        sns.violinplot(
            data=sub, x="sample", y="freq",
            order=SAMPLE_ORDER,
            palette={s: SAMPLE_PAL[s] for s in SAMPLE_ORDER},
            ax=ax,
            inner="box",
            cut=0,
            linewidth=0.7,
            density_norm="width",
        )

        n_sites = mc[mc["seq_ctx"] == ctx].groupby("sample").size()
        note = "\n".join(f"{s}: {n_sites.get(s, 0):,}" for s in SAMPLE_ORDER)
        ax.text(0.97, 0.97, note, transform=ax.transAxes,
                ha="right", va="top", fontsize=6.5, color="#555555", family="monospace")

        ax.set_title(ctx)
        ax.set_xlabel("")
        ax.set_ylim(-0.05, 1.05)
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
        if ax is axes[0]:
            ax.set_ylabel("Methylation frequency")
        else:
            ax.set_ylabel("")

    fig.suptitle("5mC Methylation Frequency by Sequence Context", y=1.02, fontsize=11)
    fig.text(0.5, -0.02, cov_note, ha="center", fontsize=7.5, color="#888888")
    plt.tight_layout()
    save(fig, "fig02_5mc_by_seq_context.png")

    # ── fig03 – chromosome-level methylation density by sequence context ──────
    # One figure per context; 10 horizontal strips, 500 kb bins, 4 sample lines.

    chrom_sizes = pd.read_csv(
        str(REF_FA) + ".fai", sep="\t", header=None,
        names=["chrom", "length", "offset", "lb", "lw"],
        usecols=["chrom", "length"],
    )
    main_chroms = (
        chrom_sizes[chrom_sizes["length"] >= 10_000_000]
        .sort_values("chrom")["chrom"].tolist()
    )
    chr_label = {c: f"Chr{i+1:02d}" for i, c in enumerate(main_chroms)}
    chr_len   = chrom_sizes.set_index("chrom")["length"]

    BIN_SIZE = 500_000
    mc_chr = mc[mc["chrom"].isin(main_chroms)].copy()
    mc_chr["bin_mb"] = (mc_chr["start"] // BIN_SIZE) * BIN_SIZE / 1e6

    binned = (
        mc_chr.groupby(["chrom", "bin_mb", "sample", "seq_ctx"])["freq"]
        .mean()
        .reset_index()
    )

    for ctx in SEQ_CTXS:
        ctx_df = binned[binned["seq_ctx"] == ctx]

        fig, axes = plt.subplots(
            len(main_chroms), 1,
            figsize=(13, 1.8 * len(main_chroms)),
            sharey=True, sharex=False,
            gridspec_kw={"hspace": 0.08},
        )

        for ax, chrom in zip(axes, main_chroms):
            clen_mb = chr_len[chrom] / 1e6
            cdf = ctx_df[ctx_df["chrom"] == chrom]

            for sample in SAMPLE_ORDER:
                sdf = cdf[cdf["sample"] == sample].sort_values("bin_mb")
                if sdf.empty:
                    continue
                ax.plot(sdf["bin_mb"], sdf["freq"],
                        color=SAMPLE_PAL[sample], lw=0.7, alpha=0.85, label=sample)

            ax.set_xlim(0, clen_mb)
            ax.set_ylim(0, 1)
            ax.set_yticks([0, 0.5, 1])
            ax.tick_params(axis="y", labelsize=6.5)
            ax.tick_params(axis="x", labelsize=7)
            ax.text(0.005, 0.90, chr_label[chrom],
                    transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top")
            if ax is not axes[-1]:
                ax.tick_params(labelbottom=False)

        axes[-1].set_xlabel("Position (Mb)")
        fig.text(-0.01, 0.5, "Mean methylation frequency",
                 va="center", rotation="vertical", fontsize=9)

        handles = [
            plt.Line2D([0], [0], color=SAMPLE_PAL[s], lw=1.5, label=s)
            for s in SAMPLE_ORDER
        ]
        axes[0].legend(handles=handles, loc="upper right",
                       fontsize=7, ncol=len(SAMPLE_ORDER), framealpha=0.7)

        fig.suptitle(f"{ctx} Methylation Density along Chromosomes", y=1.005, fontsize=11)
        fig.text(0.5, -0.01, cov_note, ha="center", fontsize=7.5, color="#888888")
        save(fig, f"fig03_chrplot_{ctx}.png")

logging.info(f"All figures saved → {OUT_DIR}/")
