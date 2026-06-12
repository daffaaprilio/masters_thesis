#!/usr/bin/env python3
"""
Variant intersection UpSet plot across four sorghum samples.

Runs bcftools isec on all four phased VCFs (no -n filter) to capture every
combination of sample overlap, parses sites.txt for the membership pattern,
and generates an UpSet plot.

Run via:
    ./docker/run.sh python3 analysis/scripts/variant_venn.py

Output: analysis/03_TAA/figures/variant_upset.png
"""

import subprocess
import tempfile
import logging
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

# ── Config ─────────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).resolve().parent.parent.parent
VCF_DIR = ROOT / "results/vcf_processing"
OUT_DIR = ROOT / "analysis/03_TAA/figures"
LOG_DIR = ROOT / "analysis/logs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]
SAMPLE_LABELS = {
    "SBC4":  "SBC4 (TAA ++, high sec.)",
    "SBC10": "SBC10 (TAA +++, low sec.)",
    "SBC11": "SBC11 (TAA −, high sec.)",
    "SBC23": "SBC23 (TAA ++, high sec.)",
}

_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(LOG_DIR / f"variant_venn_{_ts}.log"),
    ],
)

warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

plt.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "figure.dpi":     150,
})

# ── Run bcftools isec ──────────────────────────────────────────────────────────
vcfs = [str(VCF_DIR / f"{s}.phased.vcf.gz") for s in SAMPLES]
logging.info("Running bcftools isec on phased VCFs…")

memberships: list[list[str]] = []
with tempfile.TemporaryDirectory() as tmp:
    result = subprocess.run(
        ["bcftools", "isec", "-p", tmp] + vcfs,
        check=True, capture_output=True, text=True,
    )
    if result.stderr:
        logging.info(result.stderr.strip())

    sites_path = Path(tmp) / "sites.txt"
    with open(sites_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # columns: CHROM  POS  REF  ALT  <binary pattern>
            pattern = line.rstrip().split("\t")[4]
            members = [SAMPLES[i] for i, c in enumerate(pattern) if c == "1"]
            if members:
                memberships.append(members)

total = len(memberships)
logging.info(f"Total variant sites parsed: {total:,}")

# Log per-combination counts for quick review
pattern_counts = Counter(tuple(sorted(m)) for m in memberships)
for combo, cnt in sorted(pattern_counts.items(), key=lambda x: -x[1]):
    logging.info(f"  {' & '.join(combo):<45}  {cnt:>7,}")

# ── Build UpSet data ───────────────────────────────────────────────────────────
labeled = [
    [SAMPLE_LABELS[s] for s in members]
    for members in memberships
]
data = from_memberships(labeled)

# ── Plot ───────────────────────────────────────────────────────────────────────
# Let upsetplot own the figure to avoid GridSpec conflicts with a pre-created fig.
# show_counts=False avoids a matplotlib text-annotation bug (0-d array → float).
upset = UpSet(
    data,
    subset_size="count",
    sort_by="cardinality",
    show_counts=False,
    totals_plot_elements=4,
)
axes_dict = upset.plot()
fig = plt.gcf()
fig.set_size_inches(13, 6)

axes_dict["intersections"].set_title(
    f"Variant set intersections — four sorghum accessions"
    f"  (n = {total:,} sites, SNP + INDEL, phased VCFs)",
    fontsize=10, pad=8,
)

out = OUT_DIR / "variant_upset.png"
fig.savefig(out, dpi=150)
logging.info(f"Saved → {out}")
print(f"\nFigure written to: {out}")
