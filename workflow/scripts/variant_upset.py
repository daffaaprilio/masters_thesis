#!/usr/bin/env python3
"""
Variant intersection UpSet plot across four sorghum samples.

Reads the pipeline's final combined multi-sample VCF (rule annotate_all →
results/combined/all.annotated.vcf.gz), which concatenates both tracks into one
file: the SNP/indel track (SnpEff+SIFT4G, merged.sift4g.vcf.gz) and the SV track
(SnpEff, snpeff_sv/combined.annotated.vcf.gz). There is no per-group VCF to read
counts from anymore — each record's group membership is derived on the fly from
its genotypes (a sample "has" a record if its GT carries an ALT allele, i.e. HET
or HOM_ALT; "./." and "0/0" don't count), the same GT-presence semantics the
pipeline itself uses to split SV records into sample-sharing groups.

Run via:
    ./docker/run.sh python3 workflow/scripts/variant_upset.py

Output: analysis/figures/variant_upset.png
"""

import logging
import warnings
from collections import Counter
from datetime import datetime
from itertools import combinations
from pathlib import Path

import cyvcf2
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

# ── Config ─────────────────────────────────────────────────────────────────────
ROOT         = Path(__file__).resolve().parent.parent.parent
COMBINED_VCF = ROOT / "results/combined/all.annotated.vcf.gz"
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


# ── Tally intersection sizes directly from genotypes ────────────────────────────
# cyvcf2 gt_types encoding: 0=HOM_REF, 1=HET, 2=UNKNOWN (missing), 3=HOM_ALT.
HAS_ALT = (1, 3)

if not COMBINED_VCF.exists():
    raise SystemExit(
        f"Missing combined VCF: {COMBINED_VCF}\n"
        "Run the pipeline first, e.g.:\n"
        "  ./docker/run.sh snakemake results/combined/all.annotated.vcf.gz --cores <N>"
    )

logging.info(f"Reading combined VCF: {COMBINED_VCF}")
vcf = cyvcf2.VCF(str(COMBINED_VCF))
sample_idx = [vcf.samples.index(s) for s in SAMPLES]

counts: Counter = Counter()
n_snp = n_sv = n_no_member = 0

for rec in vcf:
    gt = rec.gt_types
    members = tuple(s for s, i in zip(SAMPLES, sample_idx) if gt[i] in HAS_ALT)
    if not members:
        n_no_member += 1
        continue
    counts[members] += 1
    if rec.INFO.get("SVTYPE") is not None:
        n_sv += 1
    else:
        n_snp += 1

if n_no_member:
    logging.warning(f"{n_no_member:,} record(s) had no ALT-carrying sample; skipped")

total = n_snp + n_sv
logging.info(f"Total variant records: {total:,}  ({n_snp:,} SNP/INDEL, {n_sv:,} SV)")

# ── Build UpSet data ───────────────────────────────────────────────────────────
# Every non-empty subset of SAMPLES, in canonical order — same 15 groups used
# elsewhere in the pipeline (results/variant_groups, results/sv_groups labels).
memberships: list[list[str]] = []
group_counts: list[int] = []
for k in range(1, len(SAMPLES) + 1):
    for combo in combinations(SAMPLES, k):
        n = counts.get(combo, 0)
        memberships.append([SAMPLE_LABELS[s] for s in combo])
        group_counts.append(n)
        logging.info(f"  {'_'.join(combo):<28}  {n:>10,}")

# One row per membership combination; the value is that combination's record count.
data = from_memberships(memberships, data=group_counts)

# ── Plot ───────────────────────────────────────────────────────────────────────
# Let upsetplot own the figure to avoid GridSpec conflicts with a pre-created fig.
# subset_size="sum" sums each combination's record count (one row per combination).
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
    f"  (n = {total:,} sites)",
    fontsize=10, pad=8,
)

out = OUT_DIR / "variant_upset.png"
fig.savefig(out, dpi=150)
logging.info(f"Saved → {out}")
print(f"\nFigure written to: {out}")
