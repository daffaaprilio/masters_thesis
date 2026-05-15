#!/Users/daffa/miniconda3/envs/sbi/bin/python
"""
VCF Quality Benchmark: QUAL and DP distributions.

Compares raw Clair3 VCFs (results/vcf/) vs annotated/filtered VCFs
(analysis/data/vcf/annotated/) for each sample.

Output (analysis/01_variant_landscape/figures/):
  fig0_vcf_benchmark.png  – QUAL and DP violin distributions, raw vs annotated
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────────────────────

SAMPLES  = ["SBC4", "SBC10", "SBC11", "SBC23"]
RAW_DIR  = Path(__file__).parent.parent / "results/vcf"
ANN_DIR  = Path(__file__).parent / "data/vcf/annotated"
OUT_DIR  = Path(__file__).parent / "01_variant_landscape/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAX_N = 80_000   # reservoir-sample cap per VCF

plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       9,
    "axes.labelsize":  9,
    "axes.titlesize":  10,
    "legend.fontsize": 8,
    "figure.dpi":      150,
})

COLORS = {"Raw (Clair3)": "#7fbfdf", "Annotated": "#f4a261"}

# ── Extraction ─────────────────────────────────────────────────────────────────

def _reservoir(iterable, k: int, rng: np.random.Generator):
    """Reservoir sampling (algorithm R) over an iterable."""
    reservoir = []
    for i, item in enumerate(iterable):
        if i < k:
            reservoir.append(item)
        else:
            j = int(rng.integers(0, i + 1))
            if j < k:
                reservoir[j] = item
    return reservoir


def extract_qual_dp(vcf_path: Path, max_n: int = MAX_N) -> pd.DataFrame:
    """Stream QUAL + FORMAT/DP from a VCF via bcftools; reservoir-sample to max_n rows."""
    cmd = ["bcftools", "query", "-f", r"%QUAL\t[%DP]\n", str(vcf_path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)

    rng = np.random.default_rng(42)
    sample_lines = _reservoir(proc.stdout, max_n, rng)
    proc.stdout.close()
    proc.wait()

    rows = []
    for line in sample_lines:
        parts = line.rstrip().split("\t")
        if len(parts) == 2:
            try:
                rows.append((float(parts[0]), int(parts[1])))
            except ValueError:
                pass
    return pd.DataFrame(rows, columns=["QUAL", "DP"])


# ── Load ───────────────────────────────────────────────────────────────────────

frames = []
counts: dict[str, dict[str, int]] = {}

for s in SAMPLES:
    df_raw = extract_qual_dp(RAW_DIR / f"{s}.vcf.gz")
    df_ann = extract_qual_dp(ANN_DIR / f"{s}.annotated.vcf.gz")
    df_raw["sample"] = s
    df_raw["source"] = "Raw (Clair3)"
    df_ann["sample"] = s
    df_ann["source"] = "Annotated"
    frames.extend([df_raw, df_ann])

    # True variant counts from bcftools stats
    def _count(path):
        out = subprocess.run(
            ["bcftools", "stats", str(path)],
            capture_output=True, text=True
        ).stdout
        for line in out.splitlines():
            if "number of records:" in line:
                return int(line.split()[-1])
        return 0

    counts[s] = {
        "raw": _count(RAW_DIR / f"{s}.vcf.gz"),
        "ann": _count(ANN_DIR / f"{s}.annotated.vcf.gz"),
    }
    print(f"{s}: raw={counts[s]['raw']:,}  annotated={counts[s]['ann']:,}")

df = pd.concat(frames, ignore_index=True)

# ── Figure ─────────────────────────────────────────────────────────────────────

fig, (ax_q, ax_d) = plt.subplots(2, 1, figsize=(11, 9), constrained_layout=True)
fig.suptitle("VCF Quality Benchmark: Raw Clair3 vs Annotated/Filtered",
             fontsize=12, fontweight="bold")

# ── Panel A: QUAL distribution ─────────────────────────────────────────────────

# Clip QUAL at 0.1 for log scale; cap display at 120
df_q = df.copy()
df_q["QUAL"] = df_q["QUAL"].clip(lower=0.1, upper=200)

order = SAMPLES
hue_order = ["Raw (Clair3)", "Annotated"]
palette = [COLORS[h] for h in hue_order]

sns.violinplot(
    data=df_q, x="sample", y="QUAL", hue="source",
    order=order, hue_order=hue_order, palette=palette,
    inner="quartile", linewidth=0.8, density_norm="width",
    ax=ax_q,
)
ax_q.set_yscale("log")
ax_q.set_ylim(0.1, 200)
ax_q.axhline(20, color="#e63946", ls="--", lw=1.2, label="QUAL = 20 (filter threshold)")
ax_q.set_title("A   Genotype Quality (QUAL) Distribution")
ax_q.set_xlabel("")
ax_q.set_ylabel("QUAL score (log scale)")
ax_q.xaxis.set_ticks(range(len(order)))
ax_q.set_xticklabels(order)

# Annotate variant counts above each group
for i, s in enumerate(SAMPLES):
    ax_q.text(i - 0.22, 150,
              f"n={counts[s]['raw']/1e6:.2f}M", ha="center", va="bottom",
              fontsize=7, color="#2a6f97", fontstyle="italic")
    ax_q.text(i + 0.22, 150,
              f"n={counts[s]['ann']/1e6:.2f}M", ha="center", va="bottom",
              fontsize=7, color="#c6632a", fontstyle="italic")

thresh_line = mlines.Line2D([], [], color="#e63946", ls="--", lw=1.2,
                             label="QUAL = 20 filter threshold")
handles, labels = ax_q.get_legend_handles_labels()
# Replace seaborn legend with custom one that includes threshold line
ax_q.legend(handles=handles[:2] + [thresh_line],
            labels=hue_order + ["QUAL = 20 (filter threshold)"],
            loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)

# ── Panel B: DP distribution ──────────────────────────────────────────────────

# Cap DP display at 200 so extreme outliers don't squash the violins
df_d = df.copy()
df_d["DP"] = df_d["DP"].clip(upper=200)

sns.violinplot(
    data=df_d, x="sample", y="DP", hue="source",
    order=order, hue_order=hue_order, palette=palette,
    inner="quartile", linewidth=0.8, density_norm="width",
    ax=ax_d,
)
ax_d.axhline(10,  color="#e9c46a", ls="--", lw=1.2, label="DP = 10 (min)")
ax_d.axhline(100, color="#e63946", ls="--", lw=1.2, label="DP = 100 (max)")
ax_d.set_title("B   Read Depth (DP) Distribution")
ax_d.set_xlabel("Sample")
ax_d.set_ylabel("Read depth (DP)")
ax_d.set_ylim(0, 210)
ax_d.xaxis.set_ticks(range(len(order)))
ax_d.set_xticklabels(order)

handles_d, _ = ax_d.get_legend_handles_labels()
dp_min_line = mlines.Line2D([], [], color="#e9c46a", ls="--", lw=1.2, label="DP = 10 (min filter)")
dp_max_line = mlines.Line2D([], [], color="#e63946", ls="--", lw=1.2, label="DP = 100 (max filter)")
ax_d.legend(handles=handles_d[:2] + [dp_min_line, dp_max_line],
            labels=hue_order + ["DP = 10 (min filter)", "DP = 100 (max filter)"],
            loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)

out = OUT_DIR / "fig0_vcf_benchmark.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved → {out}")
