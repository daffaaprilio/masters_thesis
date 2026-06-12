#!/usr/bin/env python3
"""
Variant landscape EDA from SnpEff per-sample stats CSVs.

Outputs (analysis/01_variant_landscape/figures/):
  fig01_variant_type.png
  fig02_effect_impact.png      – broken y-axis: 0–10% shown, 10–100% skipped
  fig03_genomic_region.png
  fig04_tstv_ratio.png
  fig05_functional_class.png
  fig06_ms_ratio.png
  fig07_indel_length.png
  fig08_chr_density.png
  fig09_base_changes_<sample>.png  (one per sample)
"""

import re
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from datetime import datetime
from io import StringIO
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).parent.parent / "data/vcf/annotated"
OUT_DIR  = Path(__file__).parent.parent / "01_variant_landscape/figures"
LOG_DIR  = Path(__file__).parent.parent / "logs"
OUT_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
_log_path  = LOG_DIR / f"variant_landscape_{_timestamp}.log"
_root = logging.getLogger()
_root.setLevel(logging.INFO)
_ch = logging.StreamHandler()
_ch.setFormatter(logging.Formatter("%(message)s"))
_root.addHandler(_ch)
_fh = logging.FileHandler(_log_path)
_fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
_root.addHandler(_fh)
logging.info(f"Log: {_log_path}")

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]

plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       9,
    "axes.labelsize":  9,
    "axes.titlesize":  10,
    "legend.fontsize": 7.5,
    "figure.dpi":      150,
})

# ── Section parsers ────────────────────────────────────────────────────────────

def split_sections(path: Path) -> dict[str, str]:
    text = path.read_text()
    matches = list(re.finditer(r'^# (.+)$', text, re.MULTILINE))
    out: dict[str, str] = {}
    for i, m in enumerate(matches):
        name  = m.group(1).strip()
        start = m.end()
        end   = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        out[name] = text[start:end].strip()
    return out


def find_section(sections: dict[str, str], keyword: str) -> str:
    kw = keyword.lower()
    for k, v in sections.items():
        if kw in k.lower():
            return v
    return ""


def parse_csv(block: str) -> pd.DataFrame:
    lines = [l for l in block.splitlines() if l.strip()]
    if not lines:
        return pd.DataFrame()
    try:
        df = pd.read_csv(StringIO("\n".join(lines)), sep=r"\s*,\s*", engine="python")
        df.columns = [c.strip() for c in df.columns]
        for col in df.select_dtypes("object").columns:
            if df[col].str.endswith("%", na=False).any():
                df[col] = df[col].str.rstrip("%").astype(float)
        return df
    except Exception:
        return pd.DataFrame()


def parse_kv(block: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",", 1)]
        if len(parts) == 2 and parts[1]:
            out[parts[0]] = parts[1].split()[0]
    return out


def parse_values_counts(block: str) -> pd.Series:
    lines = [l.strip() for l in block.splitlines() if l.strip()]
    if len(lines) < 2:
        return pd.Series(dtype=float)
    vals = [v.strip() for v in lines[0].split(",")[1:]]
    cnts = [v.strip() for v in lines[1].split(",")[1:]]
    try:
        return pd.Series({int(v): int(c) for v, c in zip(vals, cnts) if v and c})
    except Exception:
        return pd.Series(dtype=float)


def parse_base_changes(block: str) -> pd.DataFrame:
    lines = [l.strip() for l in block.splitlines() if l.strip()]
    if not lines:
        return pd.DataFrame()
    to_bases = [x.strip() for x in lines[0].split(",")[1:5]]
    rows: dict[str, dict] = {}
    for line in lines[1:]:
        parts = [x.strip() for x in line.split(",")]
        if len(parts) >= 5:
            try:
                rows[parts[0]] = {b: int(parts[i + 1]) for i, b in enumerate(to_bases)}
            except (ValueError, IndexError):
                pass
    return pd.DataFrame(rows).T if rows else pd.DataFrame()


# ── Load samples ───────────────────────────────────────────────────────────────

def load_sample(name: str) -> dict:
    path = DATA_DIR / f"{name}.stats.csv"
    secs = split_sections(path)

    fc_block = find_section(secs, "functional class")
    ms_ratio: float | None = None
    for line in fc_block.splitlines():
        if "Missense_Silent_ratio" in line:
            try:
                ms_ratio = float(line.split(",")[1].strip())
            except (IndexError, ValueError):
                pass

    return {
        "chr_rate":   parse_csv(find_section(secs, "Change rate")),
        "var_type":   parse_csv(find_section(secs, "variant")),
        "impact":     parse_csv(find_section(secs, "by impact")),
        "func_class": parse_csv(fc_block),
        "region":     parse_csv(find_section(secs, "genomic region")),
        "quality":    parse_values_counts(find_section(secs, "Quality")),
        "indel_len":  parse_values_counts(find_section(secs, "InDel")),
        "base_chg":   parse_base_changes(find_section(secs, "Base changes")),
        "tstv":       parse_kv(find_section(secs, "Ts/Tv summary")),
        "ms_ratio":   ms_ratio,
    }


data = {s: load_sample(s) for s in SAMPLES}

# ── Shared helpers ─────────────────────────────────────────────────────────────

def collect_stacked(key: str, keep: list | None = None) -> pd.DataFrame:
    frames: dict[str, pd.Series] = {}
    for s in SAMPLES:
        df = data[s][key]
        if df.empty or "Type" not in df.columns:
            continue
        ser = df.dropna(subset=["Count"]).set_index("Type")["Count"]
        # drop metadata rows (e.g. "Missense_Silent_ratio", "Change_rate");
        # use word-boundary check to avoid accidentally dropping "MODERATE"
        ser = ser[ser.index.map(
            lambda x: isinstance(x, str)
            and not re.search(r'(ratio|_rate)', x, re.IGNORECASE)
        )]
        if keep is not None:
            ser = ser.reindex(keep).fillna(0)
        frames[s] = ser
    return pd.DataFrame(frames).fillna(0).T


def pct(df: pd.DataFrame) -> pd.DataFrame:
    return df.div(df.sum(axis=1), axis=0) * 100


def stacked_bar(
    ax: plt.Axes,
    df: pd.DataFrame,
    colors: dict,
    title: str,
    ylabel: str = "%",
) -> None:
    bottom = np.zeros(len(df))
    x = np.arange(len(df))
    for col in df.columns:
        ax.bar(x, df[col], bottom=bottom, color=colors.get(col, "#aaaaaa"),
               label=col, width=0.6)
        bottom += df[col].values
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=30, ha="right")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 103)
    ax.legend(loc="upper right", framealpha=0.8)


def save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved {name}")


# ── Colour palettes ────────────────────────────────────────────────────────────

VARTYPE_ORDER  = ["SNP", "INS", "DEL"]
VARTYPE_COLORS = {"SNP": "#4c72b0", "INS": "#55a868", "DEL": "#dd8452"}

IMPACT_ORDER  = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
IMPACT_COLORS = {
    "HIGH":     "#c44e52",
    "MODERATE": "#dd8452",
    "LOW":      "#4c72b0",
    "MODIFIER": "#8172b2",
}

TOP_REGIONS = [
    "INTERGENIC", "UPSTREAM", "DOWNSTREAM", "INTRON",
    "EXON", "UTR_3_PRIME", "UTR_5_PRIME", "SPLICE_SITE_REGION",
]
REGION_COLORS = dict(zip(TOP_REGIONS, sns.color_palette("tab10", len(TOP_REGIONS))))
REGION_COLORS["OTHER"] = "#aaaaaa"

FC_ORDER  = ["MISSENSE", "SILENT", "NONSENSE"]
FC_COLORS = {"MISSENSE": "#e07b54", "SILENT": "#5b9bd5", "NONSENSE": "#c44e52"}

# ── fig01 – Variant type ───────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))
df_vt = pct(collect_stacked("var_type", keep=VARTYPE_ORDER))
df_vt = df_vt[[c for c in VARTYPE_ORDER if c in df_vt.columns]]
stacked_bar(ax, df_vt, VARTYPE_COLORS, "Variant Type Composition")
plt.tight_layout()
save(fig, "fig01_variant_type.png")

# ── fig02 – Effect impact (broken y-axis: 0–LO_MAX then HI_MIN–100) ──────────

df_imp = pct(collect_stacked("impact", keep=IMPACT_ORDER))
df_imp = df_imp[[c for c in IMPACT_ORDER if c in df_imp.columns]]

# size panels around the data: bottom shows all non-MODIFIER stacked + padding
max_small = (df_imp.drop(columns="MODIFIER", errors="ignore")
                   .sum(axis=1).max())
LO_MAX = np.ceil(max_small * 1.35)   # e.g. ~4 %
HI_MIN = 93.0                         # top panel starts here

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

# x-axis labels only on bottom panel
ax_bot.set_xticks(x)
ax_bot.set_xticklabels(df_imp.index, rotation=30, ha="right")
ax_top.tick_params(axis="x", bottom=False)

# hide the inner spines to make the break obvious
ax_top.spines.bottom.set_visible(False)
ax_bot.spines.top.set_visible(False)

# diagonal break marks at the gap
d = 0.018
kw = dict(color="k", clip_on=False, lw=1.2)
ax_top.plot([-d, +d], [-d, +d], transform=ax_top.transAxes, **kw)
ax_top.plot([1-d, 1+d], [-d, +d], transform=ax_top.transAxes, **kw)
ax_bot.plot([-d, +d], [1-d, 1+d], transform=ax_bot.transAxes, **kw)
ax_bot.plot([1-d, 1+d], [1-d, 1+d], transform=ax_bot.transAxes, **kw)

ax_bot.set_ylabel("% of effects")
ax_top.set_title("Effect Impact Distribution")
ax_top.legend(loc="lower right", framealpha=0.8, fontsize=7.5)

plt.tight_layout()

# place ⋮ on the y-axis between the two panels
bot_top = ax_bot.get_position().y1    # top edge of bottom subplot (figure coords)
top_bot = ax_top.get_position().y0   # bottom edge of top subplot
x_left  = ax_bot.get_position().x0
fig.text(x_left - 0.04, (bot_top + top_bot) / 2, "⋮",
         ha="center", va="center", fontsize=13, color="gray",
         transform=fig.transFigure, clip_on=False)

save(fig, "fig02_effect_impact.png")

# ── fig03 – Genomic region ────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))
raw_reg = collect_stacked("region")
top_cols = [c for c in TOP_REGIONS if c in raw_reg.columns]
other_cols = [c for c in raw_reg.columns if c not in TOP_REGIONS]
df_reg = raw_reg[top_cols].copy()
if other_cols:
    df_reg["OTHER"] = raw_reg[other_cols].sum(axis=1)
df_reg = pct(df_reg)
stacked_bar(ax, df_reg, REGION_COLORS, "Genomic Region Distribution")
plt.tight_layout()
save(fig, "fig03_genomic_region.png")

# ── fig04 – Ts/Tv ratio ───────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(5, 4))
tstv_vals = {s: float(data[s]["tstv"].get("Ts_Tv_ratio", 0)) for s in SAMPLES}
bars = ax.bar(range(len(SAMPLES)), list(tstv_vals.values()),
              color=sns.color_palette("Set2", len(SAMPLES)), width=0.5)
ax.set_xticks(range(len(SAMPLES)))
ax.set_xticklabels(SAMPLES, rotation=30, ha="right")
ax.axhline(2.0, color="red", ls="--", lw=1, label="Expected ≈ 2.0")
ax.set_title("Ts/Tv Ratio")
ax.set_ylabel("Ts/Tv")
ax.set_ylim(0, max(tstv_vals.values()) * 1.25)
ax.legend(fontsize=7)
for bar, val in zip(bars, tstv_vals.values()):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
            f"{val:.3f}", ha="center", va="bottom", fontsize=8)
plt.tight_layout()
save(fig, "fig04_tstv_ratio.png")

# ── fig05 – Functional class ──────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))
df_fc = pct(collect_stacked("func_class", keep=FC_ORDER))
df_fc = df_fc[[c for c in FC_ORDER if c in df_fc.columns]]
stacked_bar(ax, df_fc, FC_COLORS, "Coding Variant Functional Class")
plt.tight_layout()
save(fig, "fig05_functional_class.png")

# ── fig06 – Missense : Silent ratio ───────────────────────────────────────────

fig, ax = plt.subplots(figsize=(5, 4))
ms = {s: data[s]["ms_ratio"] for s in SAMPLES if data[s]["ms_ratio"] is not None}
bars_ms = ax.bar(range(len(ms)), list(ms.values()),
                 color=sns.color_palette("Set2", len(ms)), width=0.5)
ax.set_xticks(range(len(ms)))
ax.set_xticklabels(list(ms.keys()), rotation=30, ha="right")
ax.axhline(1.0, color="gray", ls="--", lw=1, label="ratio = 1")
ax.set_title("Missense : Silent Ratio")
ax.set_ylabel("Ratio")
ax.legend(fontsize=7)
for i, (s, v) in enumerate(ms.items()):
    ax.text(i, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
plt.tight_layout()
save(fig, "fig06_ms_ratio.png")

# ── fig07 – InDel length distribution ────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))
colors_il = sns.color_palette("Set1", len(SAMPLES))
for i, s in enumerate(SAMPLES):
    ser = data[s]["indel_len"]
    if not ser.empty:
        ser = ser[(ser.index >= 1) & (ser.index <= 30)]
        ax.plot(ser.index, ser.values, marker="o", ms=3,
                lw=1.5, label=s, color=colors_il[i], alpha=0.85)
ax.set_yscale("log")
ax.set_title("InDel Length Distribution")
ax.set_xlabel("InDel length (bp)")
ax.set_ylabel("Count (log scale)")
ax.legend(fontsize=7)
ax.set_xlim(0.5, 30.5)
plt.tight_layout()
save(fig, "fig07_indel_length.png")

# ── fig08 – Per-chromosome variant density ────────────────────────────────────

chr_rows: dict[str, pd.Series] = {}
for s in SAMPLES:
    df_c = data[s]["chr_rate"]
    if df_c.empty or "Change_rate" not in df_c.columns:
        continue
    df_c = df_c.copy()
    df_c["Chromosome"] = df_c["Chromosome"].astype(str).str.strip()
    chr_rows[s] = df_c.set_index("Chromosome")["Change_rate"].astype(float)

if chr_rows:
    heat = pd.DataFrame(chr_rows).T
    cols = sorted(heat.columns, key=lambda x: int(x) if x.isdigit() else 99)
    heat = heat[cols]
    fig, ax = plt.subplots(figsize=(9, 3.5))
    sns.heatmap(
        heat, ax=ax, cmap="YlOrRd_r",
        annot=True, fmt=".0f", linewidths=0.4,
        cbar_kws={"label": "bp / variant\n(lower = denser)"},
    )
    ax.set_title("Per-Chromosome Variant Density")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    plt.tight_layout()
    save(fig, "fig08_chr_density.png")

# ── fig09–12 – Base substitution spectra (one per sample) ────────────────────

for idx, s in enumerate(SAMPLES, start=9):
    fig, ax = plt.subplots(figsize=(5, 4))
    df_bc = data[s]["base_chg"].copy()
    if df_bc.empty:
        ax.axis("off")
        save(fig, f"fig{idx:02d}_base_changes_{s}.png")
        continue
    df_bc = df_bc.astype(float)
    for b in df_bc.index:
        if b in df_bc.columns:
            df_bc.loc[b, b] = np.nan
    annot = df_bc.map(lambda x: "" if pd.isna(x) else f"{x:.0f}")
    sns.heatmap(
        df_bc, ax=ax, cmap="Blues", annot=annot, fmt="",
        linewidths=0.5, linecolor="white",
        cbar_kws={"shrink": 0.8, "label": "count"},
        mask=df_bc.isna(),
    )
    ax.set_title(f"Base Substitution Spectrum – {s}")
    ax.set_xlabel("To base")
    ax.set_ylabel("From base")
    plt.tight_layout()
    save(fig, f"fig{idx:02d}_base_changes_{s}.png")

logging.info(f"All figures saved -> {OUT_DIR}/")
