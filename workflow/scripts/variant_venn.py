#!/usr/bin/env python3
"""
4-way Venn diagram of variant site intersections across four sorghum accessions.

Runs bcftools isec on the four phased VCFs (no -n filter), builds per-sample
sets from sites.txt, and draws a 4-way Venn diagram using the venn package.

Run via:
    ./docker/run.sh python3 analysis/scripts/variant_venn.py

Output: analysis/03_TAA/figures/variant_venn.png
"""

import subprocess
import tempfile
import logging
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from venn import venn

# ── Config ─────────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).resolve().parent.parent.parent
VCF_DIR = ROOT / "results/vcf_processing"
OUT_DIR = ROOT / "analysis/03_TAA/figures"
LOG_DIR = ROOT / "analysis/logs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]
SAMPLE_LABELS = {
    "SBC4":  "SBC4 (TAA ++)",
    "SBC10": "SBC10 (TAA +++)",
    "SBC11": "SBC11 (TAA −)",
    "SBC23": "SBC23 (TAA ++)",
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

plt.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      9,
    "axes.titlesize": 10,
    "figure.dpi":     150,
})

# ── Run bcftools isec ──────────────────────────────────────────────────────────
vcfs = [str(VCF_DIR / f"{s}.phased.vcf.gz") for s in SAMPLES]
logging.info("Running bcftools isec on phased VCFs…")

site_sets: dict[str, set] = {s: set() for s in SAMPLES}
with tempfile.TemporaryDirectory() as tmp:
    result = subprocess.run(
        ["bcftools", "isec", "-p", tmp] + vcfs,
        check=True, capture_output=True, text=True,
    )
    if result.stderr:
        logging.info(result.stderr.strip())

    with open(Path(tmp) / "sites.txt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # columns: CHROM  POS  REF  ALT  <binary pattern>
            parts   = line.rstrip().split("\t")
            site_id = (parts[0], parts[1], parts[2], parts[3])
            for i, c in enumerate(parts[4]):
                if c == "1":
                    site_sets[SAMPLES[i]].add(site_id)

total = len(set().union(*site_sets.values()))
logging.info(f"Total unique variant sites: {total:,}")
for s, st in site_sets.items():
    logging.info(f"  {s}: {len(st):,} sites")

# ── Plot ───────────────────────────────────────────────────────────────────────
labeled = {SAMPLE_LABELS[s]: site_sets[s] for s in SAMPLES}

fig, ax = plt.subplots(figsize=(10, 8))
venn(labeled, ax=ax, fontsize=10, cmap="tab10")

ax.set_title(
    f"Variant site intersections — four sorghum accessions\n"
    f"n = {total:,} unique sites  (SNP + INDEL, phased VCFs)",
    fontsize=10, pad=12,
)

out = OUT_DIR / "variant_venn.png"
fig.savefig(out, bbox_inches="tight", dpi=150)
logging.info(f"Saved → {out}")
print(f"\nFigure written to: {out}")
