#!/usr/bin/env python3
"""
Variant intersection UpSet plot across four sorghum samples.

Reads the intersected variant groups produced by the main pipeline
(rule intersect_group → results/variant_groups/{group}.vcf.gz). Each group VCF
is one UpSet cell: its member samples define the membership and its record count
(read straight from the VCF index — no decompression) is the intersection size.
This script does NOT run bcftools isec itself; run the pipeline first so the
group VCFs exist.

Run via:
    ./docker/run.sh python3 workflow/scripts/variant_upset.py

Output: analysis/figures/variant_upset.png
"""

import subprocess
import logging
import warnings
from datetime import datetime
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

# ── Config ─────────────────────────────────────────────────────────────────────
ROOT         = Path(__file__).resolve().parent.parent.parent
VARGROUP_DIR = ROOT / "results/variant_groups"
OUT_DIR      = ROOT / "analysis/figures"
LOG_DIR      = ROOT / "analysis/logs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Canonical sample order — must match the pipeline (workflow/Snakefile SAMPLES),
# because group labels are the member samples joined with "_" in this order.
SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]
SAMPLE_LABELS = {
    "SBC4":  "SBC4 (TAA High, Juice ++)",
    "SBC10": "SBC10 (TAA Low, Juice +++)",
    "SBC11": "SBC11 (TAA High, Juice -)",
    "SBC23": "SBC23 (TAA High, Juice ++)",
}

_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(LOG_DIR / f"variant_upset_{_ts}.log"),
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


# ── Read intersection sizes from the pipeline's group VCFs ──────────────────────
def count_records(vcf: Path) -> int:
    """Number of variant records in a bgzipped+indexed VCF, read from its index."""
    result = subprocess.run(
        ["bcftools", "index", "--nrecords", str(vcf)],
        check=True, capture_output=True, text=True,
    )
    return int(result.stdout.strip())


logging.info(f"Reading intersected variant groups from {VARGROUP_DIR}")

memberships: list[list[str]] = []
counts: list[int] = []
missing: list[Path] = []

# Every non-empty subset of SAMPLES, in canonical order — mirrors VARGROUPS in
# workflow/Snakefile so the {group}.vcf.gz filenames line up.
for k in range(1, len(SAMPLES) + 1):
    for combo in combinations(range(len(SAMPLES)), k):
        members = [SAMPLES[i] for i in combo]
        label   = "_".join(members)
        vcf     = VARGROUP_DIR / f"{label}.vcf.gz"
        if not vcf.exists():
            missing.append(vcf)
            continue
        n = count_records(vcf)
        memberships.append([SAMPLE_LABELS[s] for s in members])
        counts.append(n)
        logging.info(f"  {label:<28}  {n:>10,}")

if missing:
    listing = "\n".join(f"  - {p}" for p in missing)
    raise SystemExit(
        f"Missing {len(missing)} group VCF(s):\n{listing}\n"
        "Run the pipeline first to produce them, e.g.:\n"
        "  ./docker/run.sh snakemake <results/variant_groups/{group}.vcf.gz ...>"
    )

total = sum(counts)
logging.info(f"Total variant sites across all groups: {total:,}")

# ── Build UpSet data ───────────────────────────────────────────────────────────
# One row per membership combination; the value is that combination's site count.
data = from_memberships(memberships, data=counts)

# ── Plot ───────────────────────────────────────────────────────────────────────
# Let upsetplot own the figure to avoid GridSpec conflicts with a pre-created fig.
# subset_size="sum" sums each combination's site count (one row per combination).
upset = UpSet(
    data,
    subset_size="sum",
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
