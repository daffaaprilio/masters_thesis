#!/Users/daffa/miniconda3/envs/sbi/bin/python
"""
Variant landscape EDA from SnpEff per-sample stats CSVs.

Outputs (analysis/figures/):
  fig1_variant_overview.png   – variant type, impact, genomic region, Ts/Tv
  fig2_mutation_patterns.png  – functional class, MS ratio, InDel len, chr density
  fig3_base_changes.png       – per-sample base substitution spectra
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).parent / "data/vcf/annotated"
OUT_DIR  = Path(__file__).parent / "figures"
OUT_DIR.mkdir(exist_ok=True)

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
    """Parse comma-separated block; strips spaces around commas and % from values."""
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
    """Build samples×types count DataFrame. Drops ratio/rate metadata rows."""
    frames: dict[str, pd.Series] = {}
    for s in SAMPLES:
        df = data[s][key]
        if df.empty or "Type" not in df.columns:
            continue
        ser = df.dropna(subset=["Count"]).set_index("Type")["Count"]
        ser = ser[ser.index.map(
            lambda x: isinstance(x, str)
            and not any(w in x.lower() for w in ("ratio", "rate"))
        )]
        if keep is not None:
            ser = ser.reindex(keep).fillna(0)
        frames[s] = ser
    return pd.DataFrame(frames).fillna(0).T   # samples × types


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


# ── Figure 1: Variant overview (2×2) ──────────────────────────────────────────

VARTYPE_ORDER  = ["SNP", "INS", "DEL"]
VARTYPE_COLORS = {"SNP": "#4c72b0", "INS": "#55a868", "DEL": "#dd8452"}

IMPACT_ORDER   = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
IMPACT_COLORS  = {
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

fig1, axes = plt.subplots(2, 2, figsize=(12, 8))
fig1.suptitle("Variant Landscape Overview", fontsize=13, fontweight="bold")

# A – variant type
df_vt = pct(collect_stacked("var_type", keep=VARTYPE_ORDER))
df_vt = df_vt[[c for c in VARTYPE_ORDER if c in df_vt.columns]]
stacked_bar(axes[0, 0], df_vt, VARTYPE_COLORS, "A   Variant Type Composition")

# B – effect impact
df_imp = pct(collect_stacked("impact", keep=IMPACT_ORDER))
df_imp = df_imp[[c for c in IMPACT_ORDER if c in df_imp.columns]]
stacked_bar(axes[0, 1], df_imp, IMPACT_COLORS, "B   Effect Impact Distribution")

# C – genomic region (group minor categories as OTHER)
raw_reg = collect_stacked("region")
top_cols = [c for c in TOP_REGIONS if c in raw_reg.columns]
other_cols = [c for c in raw_reg.columns if c not in TOP_REGIONS]
df_reg = raw_reg[top_cols].copy()
if other_cols:
    df_reg["OTHER"] = raw_reg[other_cols].sum(axis=1)
df_reg = pct(df_reg)
stacked_bar(axes[1, 0], df_reg, REGION_COLORS, "C   Genomic Region Distribution")

# D – Ts/Tv ratio
tstv_vals = {s: float(data[s]["tstv"].get("Ts_Tv_ratio", 0)) for s in SAMPLES}
ax = axes[1, 1]
bars = ax.bar(range(len(SAMPLES)), list(tstv_vals.values()),
              color=sns.color_palette("Set2", len(SAMPLES)), width=0.5)
ax.set_xticks(range(len(SAMPLES)))
ax.set_xticklabels(SAMPLES, rotation=30, ha="right")
ax.axhline(2.0, color="red", ls="--", lw=1, label="Expected ≈ 2.0")
ax.set_title("D   Ts/Tv Ratio")
ax.set_ylabel("Ts/Tv")
ax.set_ylim(0, max(tstv_vals.values()) * 1.25)
ax.legend(fontsize=7)
for bar, val in zip(bars, tstv_vals.values()):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
            f"{val:.3f}", ha="center", va="bottom", fontsize=8)

plt.tight_layout()
fig1.savefig(OUT_DIR / "fig1_variant_overview.png", dpi=150, bbox_inches="tight")
plt.close(fig1)
print("fig1_variant_overview.png")

# ── Figure 2: Mutation patterns (2×2) ─────────────────────────────────────────

FC_ORDER  = ["MISSENSE", "SILENT", "NONSENSE"]
FC_COLORS = {"MISSENSE": "#e07b54", "SILENT": "#5b9bd5", "NONSENSE": "#c44e52"}

fig2, axes2 = plt.subplots(2, 2, figsize=(12, 9))
fig2.suptitle("Variant Mutation Patterns", fontsize=13, fontweight="bold")

# A – functional class (coding variants only)
df_fc = pct(collect_stacked("func_class", keep=FC_ORDER))
df_fc = df_fc[[c for c in FC_ORDER if c in df_fc.columns]]
stacked_bar(axes2[0, 0], df_fc, FC_COLORS, "A   Coding Variant Functional Class")

# B – Missense : Silent ratio
ms = {s: data[s]["ms_ratio"] for s in SAMPLES if data[s]["ms_ratio"] is not None}
ax_ms = axes2[0, 1]
bars_ms = ax_ms.bar(range(len(ms)), list(ms.values()),
                    color=sns.color_palette("Set2", len(ms)), width=0.5)
ax_ms.set_xticks(range(len(ms)))
ax_ms.set_xticklabels(list(ms.keys()), rotation=30, ha="right")
ax_ms.axhline(1.0, color="gray", ls="--", lw=1, label="ratio = 1")
ax_ms.set_title("B   Missense : Silent Ratio")
ax_ms.set_ylabel("Ratio")
ax_ms.legend(fontsize=7)
for i, (s, v) in enumerate(ms.items()):
    ax_ms.text(i, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=8)

# C – InDel length distribution (log scale, 1–30 bp)
ax_il = axes2[1, 0]
colors_il = sns.color_palette("Set1", len(SAMPLES))
for i, s in enumerate(SAMPLES):
    ser = data[s]["indel_len"]
    if not ser.empty:
        ser = ser[(ser.index >= 1) & (ser.index <= 30)]
        ax_il.plot(ser.index, ser.values, marker="o", ms=3,
                   lw=1.5, label=s, color=colors_il[i], alpha=0.85)
ax_il.set_yscale("log")
ax_il.set_title("C   InDel Length Distribution")
ax_il.set_xlabel("InDel length (bp)")
ax_il.set_ylabel("Count (log scale)")
ax_il.legend(fontsize=7)
ax_il.set_xlim(0.5, 30.5)

# D – per-chromosome variant density heatmap
chr_rows: dict[str, pd.Series] = {}
for s in SAMPLES:
    df_c = data[s]["chr_rate"]
    if df_c.empty or "Change_rate" not in df_c.columns:
        continue
    df_c = df_c.copy()
    df_c["Chromosome"] = df_c["Chromosome"].astype(str).str.strip()
    chr_rows[s] = df_c.set_index("Chromosome")["Change_rate"].astype(float)

if chr_rows:
    heat = pd.DataFrame(chr_rows).T   # samples × chromosomes
    cols = sorted(heat.columns, key=lambda x: int(x) if x.isdigit() else 99)
    heat = heat[cols]
    sns.heatmap(
        heat, ax=axes2[1, 1], cmap="YlOrRd_r",
        annot=True, fmt=".0f", linewidths=0.4,
        cbar_kws={"label": "bp / variant\n(lower = denser)"},
    )
    axes2[1, 1].set_title("D   Per-Chromosome Variant Density")
    axes2[1, 1].set_xlabel("Chromosome")
    axes2[1, 1].set_ylabel("")
    axes2[1, 1].tick_params(axis="x", rotation=0)

plt.tight_layout()
fig2.savefig(OUT_DIR / "fig2_mutation_patterns.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
print("fig2_mutation_patterns.png")

# ── Figure 3: Base substitution spectra (2×2, one per sample) ─────────────────

fig3, axes3 = plt.subplots(2, 2, figsize=(10, 8))
fig3.suptitle("Base Substitution Spectra (SNP counts)", fontsize=13, fontweight="bold")

for idx, s in enumerate(SAMPLES):
    ax = axes3.flatten()[idx]
    df_bc = data[s]["base_chg"].copy()
    if df_bc.empty:
        ax.axis("off")
        continue

    df_bc = df_bc.astype(float)
    for b in df_bc.index:
        if b in df_bc.columns:
            df_bc.loc[b, b] = np.nan   # mask self-substitution diagonal

    annot = df_bc.map(lambda x: "" if pd.isna(x) else f"{x:.0f}")
    sns.heatmap(
        df_bc, ax=ax, cmap="Blues", annot=annot, fmt="",
        linewidths=0.5, linecolor="white",
        cbar_kws={"shrink": 0.8, "label": "count"},
        mask=df_bc.isna(),
    )
    ax.set_title(s)
    ax.set_xlabel("To base")
    ax.set_ylabel("From base")

plt.tight_layout()
fig3.savefig(OUT_DIR / "fig3_base_changes.png", dpi=150, bbox_inches="tight")
plt.close(fig3)
print("fig3_base_changes.png")

print(f"\nAll figures saved -> {OUT_DIR}/")
