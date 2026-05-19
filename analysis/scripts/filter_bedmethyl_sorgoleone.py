#!/usr/bin/env python3
"""
Filter bedMethyl files to regions corresponding to sorgoleone
biosynthesis pathway genes and their co-expressed homologs.

Regions are taken from the GFF3 annotation by gene ID, with an optional
upstream/downstream flank. bedtools intersect is used for efficient
filtering of the large bedMethyl files.

Usage:
    python filter_bedmethyl_sorgoleone.py \
        --gff resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff \
        --bedmethyl resources/bedmethyl/SBC4.filtered.bed resources/bedmethyl/SBC10.filtered.bed ... \
        --outdir analysis/data/bedmethyl_sorgoleone \
        [--flank 2000]
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Sorgoleone pathway genes + co-expressed homologs
# Source: analysis/04_sorgoleone/sorgoleone.md
# --------------------------------------------------------------------------- #

# Gene IDs (numeric) and their labels
PATHWAY_GENES = {
    8066368: "SbDES2",
    8079957: "SbDES3",
    8080259: "SbOMT3",
    8081692: "SbCYP71AM1",
}

HOMOLOGS = {
    110435045: "SbDES2-homolog",    # SORBI_3004G260800, blastp 89.38
    8072903:   "SbDES3-homolog-1",  # SORBI_3008G002800, blastp 90.81
    8055482:   "SbDES3-homolog-2",  # SORBI_3008G003200, blastp 84.22
    8079958:   "SbDES3-homolog-3",  # SORBI_3005G002800, blastp 86.56
    8076922:   "SbOMT3-homolog-1",  # SORBI_3005G086600, blastp 92.51
    110436225: "SbOMT3-homolog-2",  # SORBI_3006G008000, blastp 97.86
    8085153:   "SbOMT3-homolog-3",  # SORBI_3007G074800, blastp 70
}

ALL_GENES = {**PATHWAY_GENES, **HOMOLOGS}
LOC_NAMES = {f"LOC{gid}": label for gid, label in ALL_GENES.items()}

GFF3_COLS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]


def parse_gff_attribute(attrs, key):
    for field in attrs.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def extract_gene_regions(gff_path, loc_names, flank=0):
    """Return a BED-format DataFrame for the requested LOC gene names."""
    print(f"Reading GFF3: {gff_path}")
    gff = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=GFF3_COLS)

    genes = gff[gff["feature"] == "gene"].copy()
    genes["loc_name"] = genes["attributes"].apply(lambda a: parse_gff_attribute(a, "Name"))

    selected = genes[genes["loc_name"].isin(loc_names)].copy()
    selected["label"] = selected["loc_name"].map(loc_names)

    missing = set(loc_names) - set(selected["loc_name"])
    if missing:
        print(f"  WARNING: {len(missing)} gene(s) not found in GFF3: {', '.join(sorted(missing))}")

    def apply_flank(row, flank):
        start = row["start"] - 1  # GFF3 → BED 0-based
        end   = row["end"]
        if row["strand"] == "+":
            bed_start = max(0, start - flank)   # upstream = lower coord
            bed_end   = end                      # no downstream flank
        else:
            bed_start = start                    # no downstream flank
            bed_end   = end + flank              # upstream = higher coord
        return pd.Series({"bed_start": int(bed_start), "bed_end": int(bed_end)})

    selected[["bed_start", "bed_end"]] = selected.apply(
        apply_flank, flank=flank, axis=1
    )
    
    return selected[["seqname", "bed_start", "bed_end", "label", "loc_name"]]


def filter_bedmethyl(bedmethyl_path, regions_bed_path, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["bedtools", "intersect", "-a", str(bedmethyl_path), "-b", str(regions_bed_path), "-u"]
    with open(output_path, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True, help="Path to GFF3 annotation file")
    parser.add_argument("--bedmethyl", nargs="+", required=True, help="bedMethyl file(s) to filter")
    parser.add_argument("--outdir", required=True, help="Output directory for filtered files")
    parser.add_argument("--flank", type=int, default=2000,
                        help="Upstream/downstream flank added to each gene region in bp (default: 2000)")
    args = parser.parse_args()

    regions = extract_gene_regions(args.gff, LOC_NAMES, flank=args.flank)

    print(f"\nFound {len(regions)}/{len(LOC_NAMES)} gene regions (flank={args.flank} bp):")
    print(regions.to_string(index=False))
    print()

    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tmp:
        regions[["seqname", "bed_start", "bed_end"]].to_csv(tmp, sep="\t", header=False, index=False)
        regions_bed_path = Path(tmp.name)

    outdir = Path(args.outdir)
    for raw_path in args.bedmethyl:
        path = Path(raw_path)
        # Derive a clean sample name: SBC4.filtered.bed -> SBC4
        sample = path.name.split(".")[0]
        output_path = outdir / f"{sample}.sorgoleone.bed"
        print(f"Filtering {path.name} -> {output_path.name} ...")
        filter_bedmethyl(path, regions_bed_path, output_path)
        n_lines = int(subprocess.check_output(["wc", "-l", str(output_path)]).split()[0])
        print(f"  {n_lines:,} methylation sites retained")

    regions_bed_path.unlink()
    print("\nDone.")


if __name__ == "__main__":
    main()