#!/usr/bin/env python3
"""
Annotate DMRs with overlapping genomic features.

Builds per-feature BED files from the GFF3:
  promoter  — strand-aware upstream flank before TSS
  exon      — exon features from GFF3
  CDS       — CDS features from GFF3
  intron    — gene body minus exon union (bedtools subtract)

Then intersects each DMR against all feature BEDs using bedtools and
assembles an annotated output table.

By default the GFF3 is scoped to the 11 sorgoleone target genes.
Pass --all-genes to annotate against every gene in the GFF3 (e.g. for TAA).

Usage — sorgoleone (default):
    python analysis/scripts/annotate_DMR.py \\
        --dmr    analysis/data/sorgoleone_DMR/DMR_all_pairs_combined.tsv \\
        --gff    resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff \\
        --outdir analysis/data/sorgoleone_DMR \\
        [--flank 2000]

Usage — TAA (all genes):
    python analysis/scripts/annotate_DMR.py \\
        --dmr       analysis/data/taa_DMR/DMR_all_combined.tsv \\
        --gff       resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff \\
        --outdir    analysis/data/taa_DMR \\
        --all-genes \\
        [--flank 2000]
"""

import argparse
import logging
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Target genes (same as filter_bedmethyl_sorgoleone.py)
# --------------------------------------------------------------------------- #

PATHWAY_GENES = {
    8066368: "SbDES2",
    8079957: "SbDES3",
    8080259: "SbOMT3",
    8081692: "SbCYP71AM1",
}
HOMOLOGS = {
    110435045: "SbDES2-homolog",
    8072903:   "SbDES3-homolog-1",
    8055482:   "SbDES3-homolog-2",
    8079958:   "SbDES3-homolog-3",
    8076922:   "SbOMT3-homolog-1",
    110436225: "SbOMT3-homolog-2",
    8085153:   "SbOMT3-homolog-3",
}
LOC_NAMES    = {f"LOC{gid}": label for gid, label in {**PATHWAY_GENES, **HOMOLOGS}.items()}
LABEL_TO_EGI = {label: str(gid) for gid, label in {**PATHWAY_GENES, **HOMOLOGS}.items()}

GFF3_COLS = ["seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes"]

FEATURE_ORDER = ["promoter", "exon", "CDS", "intron"]


# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #

def setup_logging(log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter("%(message)s"))
    fh = logging.FileHandler(log_path)
    fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                                      datefmt="%Y-%m-%d %H:%M:%S"))
    root.addHandler(ch)
    root.addHandler(fh)


# --------------------------------------------------------------------------- #
# GFF3 parsing
# --------------------------------------------------------------------------- #

def parse_attr(attrs, key):
    for field in attrs.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def build_feature_beds(gff_path, loc_names, flank, feature_dir):
    """
    Parse GFF3 and write one BED file per feature type.

    When loc_names is None all genes in the GFF3 are used; otherwise only genes
    whose Name attribute appears in loc_names are included.

    Parent-child traversal:
      gene (Name=LOC...) → mRNA (Parent=gene-LOC...) → exon / CDS (Parent=mRNA-id)

    Introns are computed per gene with:
      bedtools subtract -a <gene_body.bed> -b <exons_merged.bed>

    Returns dict {feature_type: Path}.
    """
    logging.info(f"Parsing GFF3: {gff_path}")
    gff = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=GFF3_COLS)
    gff["gff_id"]     = gff["attributes"].apply(lambda a: parse_attr(a, "ID"))
    gff["gff_parent"] = gff["attributes"].apply(lambda a: parse_attr(a, "Parent"))

    # ── Step 1: target gene rows ──────────────────────────────────────────
    gene_rows = gff[gff["feature"] == "gene"].copy()
    gene_rows["loc_name"] = gene_rows["attributes"].apply(lambda a: parse_attr(a, "Name"))

    if loc_names is None:
        # All-genes mode: use every gene; label = Name attribute (LOC... or raw ID)
        target = gene_rows.copy()
        target["label"] = target["loc_name"].where(
            target["loc_name"].notna(), target["gff_id"])
    else:
        target = gene_rows[gene_rows["loc_name"].isin(loc_names)].copy()
        target["label"] = target["loc_name"].map(loc_names)

        if len(target) < len(loc_names):
            missing = set(loc_names) - set(target["loc_name"])
            logging.warning(f"{len(missing)} gene(s) not found: {', '.join(sorted(missing))}")

    gene_id_to_label = dict(zip(target["gff_id"], target["label"]))
    logging.info(f"  {len(target)} target genes found\n")

    # ── Step 2: mRNA children ─────────────────────────────────────────────
    mrna_rows = gff[
        (gff["feature"] == "mRNA") &
        (gff["gff_parent"].isin(gene_id_to_label))
    ].copy()
    mrna_id_to_label = dict(zip(mrna_rows["gff_id"],
                                mrna_rows["gff_parent"].map(gene_id_to_label)))

    # ── Step 3: exon / CDS children of target mRNAs ──────────────────────
    child_rows = gff[gff["gff_parent"].isin(mrna_id_to_label)].copy()
    child_rows["label"] = child_rows["gff_parent"].map(mrna_id_to_label)
    # Convert GFF3 1-based → BED 0-based
    child_rows["bed_start"] = child_rows["start"] - 1
    child_rows["bed_end"]   = child_rows["end"]

    feature_dir = Path(feature_dir)
    feature_dir.mkdir(parents=True, exist_ok=True)
    beds = {}

    # ── Promoter ──────────────────────────────────────────────────────────
    prom_rows = []
    for _, g in target.iterrows():
        gs = g["start"] - 1   # 0-based
        ge = g["end"]
        if g["strand"] == "+":
            ps, pe = max(0, gs - flank), gs
        else:
            ps, pe = ge, ge + flank
        if pe > ps:
            prom_rows.append({"seqname": g["seqname"], "start": ps,
                               "end": pe, "label": g["label"]})

    prom_df  = pd.DataFrame(prom_rows)
    prom_path = feature_dir / "promoter.bed"
    prom_df.sort_values(["seqname", "start"]).to_csv(
        prom_path, sep="\t", header=False, index=False)
    beds["promoter"] = prom_path
    logging.info(f"  promoter : {len(prom_df):>4} regions written → {prom_path.name}")

    # ── Exon / CDS ────────────────────────────────────────────────────────
    for feat in ("exon", "CDS"):
        rows = child_rows[child_rows["feature"] == feat]
        if rows.empty:
            logging.warning(f"  {feat}: no features found in GFF3 for target genes")
            continue
        bed = rows[["seqname", "bed_start", "bed_end", "label"]].drop_duplicates()
        bed = bed.sort_values(["seqname", "bed_start"])
        path = feature_dir / f"{feat}.bed"
        bed.to_csv(path, sep="\t", header=False, index=False)
        beds[feat] = path
        logging.info(f"  {feat:<9}: {len(bed):>4} regions written → {path.name}")

    # ── Intron: gene body − merged exons, per gene ────────────────────────
    exon_rows = child_rows[child_rows["feature"] == "exon"]
    intron_records = []

    if loc_names is None:
        logging.info("  intron   : skipped in --all-genes mode (would require ~34k subprocesses)")
        return beds

    for _, g in target.iterrows():
        gid    = g["gff_id"]
        label  = g["label"]
        gs0    = g["start"] - 1
        ge     = g["end"]

        # Exons for this gene (all mRNA isoforms)
        gene_mrna_ids = set(mrna_rows[mrna_rows["gff_parent"] == gid]["gff_id"])
        g_exons = exon_rows[exon_rows["gff_parent"].isin(gene_mrna_ids)]

        if g_exons.empty:
            logging.warning(f"  {label}: no exons found, skipping intron calculation")
            continue

        # Write gene body and exon BEDs to temp files, then bedtools subtract
        with tempfile.NamedTemporaryFile(suffix="_gbody.bed", mode="w",
                                         delete=False) as f_gene, \
             tempfile.NamedTemporaryFile(suffix="_exons.bed", mode="w",
                                         delete=False) as f_exon:
            gene_tmp = Path(f_gene.name)
            exon_tmp = Path(f_exon.name)
            f_gene.write(f"{g['seqname']}\t{gs0}\t{ge}\n")
            for _, ex in g_exons.iterrows():
                f_exon.write(f"{ex['seqname']}\t{ex['bed_start']}\t{ex['bed_end']}\n")

        result = subprocess.run(
            ["bedtools", "subtract", "-a", str(gene_tmp), "-b", str(exon_tmp)],
            capture_output=True, text=True, check=True,
        )
        gene_tmp.unlink()
        exon_tmp.unlink()

        for line in result.stdout.splitlines():
            if line:
                chrom, s, e = line.split("\t")
                intron_records.append({
                    "seqname": chrom, "start": int(s), "end": int(e), "label": label,
                })

    if intron_records:
        intron_df  = pd.DataFrame(intron_records).sort_values(["seqname", "start"])
        intron_path = feature_dir / "intron.bed"
        intron_df.to_csv(intron_path, sep="\t", header=False, index=False)
        beds["intron"] = intron_path
        logging.info(f"  intron   : {len(intron_df):>4} regions written → {intron_path.name}")
    else:
        logging.warning("  intron: no intron regions computed")

    return beds


# --------------------------------------------------------------------------- #
# DMR annotation
# --------------------------------------------------------------------------- #

def intersect_feature(dmr_bed_path, feature_bed_path, feature_type):
    """
    Run bedtools intersect -wa -wb and return a dict mapping
    (chr, start, end) → set of gene labels that overlap in this feature.
    """
    result = subprocess.run(
        ["bedtools", "intersect",
         "-a", str(dmr_bed_path),
         "-b", str(feature_bed_path),
         "-wa", "-wb"],
        capture_output=True, text=True, check=True,
    )
    hits = {}
    for line in result.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        # DMR cols: 0=chr 1=start 2=end  |  feature cols: 3=chr 4=start 5=end 6=label
        key   = (parts[0], int(parts[1]), int(parts[2]))
        label = parts[6] if len(parts) > 6 else "unknown"
        hits.setdefault(key, set()).add(label)

    n = sum(len(v) for v in hits.values())
    logging.info(f"  {feature_type:<10}: {len(hits):>3} DMRs overlapping "
                 f"({n} total feature hits)")
    return hits


def annotate_dmrs(dmr_df, feature_beds, outdir, label_to_egi=None):
    """
    Intersect DMRs against each feature BED and attach annotation columns.

    Added columns:
      features   — semicolon-separated overlapping feature types, or 'intergenic'
      gene_label — semicolon-separated gene name(s)
      egi        — semicolon-separated gene IDs (only when label_to_egi is provided)
    """
    # Write DMR BED (no header)
    dmr_bed_path = Path(outdir) / "_dmrs_tmp.bed"
    dmr_df[["chr", "start", "end"]].to_csv(
        dmr_bed_path, sep="\t", header=False, index=False)

    # Collect overlaps per feature type
    all_hits = {}   # feature_type → {key → set(labels)}
    for feat in FEATURE_ORDER:
        if feat not in feature_beds:
            continue
        all_hits[feat] = intersect_feature(dmr_bed_path, feature_beds[feat], feat)

    dmr_bed_path.unlink()

    # Build annotation per DMR row
    features_col   = []
    gene_label_col = []
    egi_col        = [] if label_to_egi is not None else None

    for _, row in dmr_df.iterrows():
        key = (row["chr"], int(row["start"]), int(row["end"]))
        feat_list  = []
        label_set  = set()

        for feat in FEATURE_ORDER:
            if feat in all_hits and key in all_hits[feat]:
                feat_list.append(feat)
                label_set |= all_hits[feat][key]

        labels = sorted(label_set)
        features_col.append(";".join(feat_list) if feat_list else "intergenic")
        gene_label_col.append(";".join(labels))
        if egi_col is not None:
            egi_col.append(";".join(label_to_egi.get(lbl, "") for lbl in labels))

    dmr_df = dmr_df.copy()
    dmr_df["features"]   = features_col
    dmr_df["gene_label"] = gene_label_col
    if egi_col is not None:
        dmr_df["egi"] = egi_col
    return dmr_df


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #

def print_summary(dmr_df):
    logging.info("\n── Feature overlap summary ──────────────────────────────────")
    for feat in FEATURE_ORDER + ["intergenic"]:
        mask = dmr_df["features"].str.contains(feat, regex=False)
        logging.info(f"  {feat:<12}: {mask.sum():>3} DMRs")

    logging.info("\n── Per-gene DMR count ───────────────────────────────────────")
    gene_counts = (
        dmr_df["gene_label"]
        .str.split(";").explode()
        .replace("", pd.NA).dropna()
        .value_counts()
        .sort_index()
    )
    for gene, n in gene_counts.items():
        logging.info(f"  {gene:<24}: {n:>3} DMRs")

    logging.info("\n── Per-pair DMR count ───────────────────────────────────────")
    pair_counts = dmr_df.groupby(["sample_a", "sample_b"]).size()
    for (sa, sb), n in pair_counts.items():
        logging.info(f"  {sa} vs {sb:<8}: {n:>3} DMRs")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dmr", default="analysis/data/sorgoleone_DMR/DMR_all_pairs_combined.tsv",
                        help="Combined DMR TSV from run_dss_dmr.R")
    parser.add_argument("--gff", default="resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff",
                        help="GFF3 annotation file")
    parser.add_argument("--outdir", default="analysis/data/sorgoleone_DMR",
                        help="Output directory (annotated TSV written here)")
    parser.add_argument("--out", default=None,
                        help="Output TSV filename (default: DMR_annotated.tsv inside --outdir)")
    parser.add_argument("--flank", type=int, default=2000,
                        help="Promoter flank in bp upstream of TSS (default: 2000)")
    parser.add_argument("--all-genes", action="store_true",
                        help="Annotate against all genes in the GFF3 instead of sorgoleone targets")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = Path(__file__).parent.parent / "logs"
    setup_logging(log_dir / f"annotate_DMR_{timestamp}.log")

    logging.info(f"DMR input : {args.dmr}")
    logging.info(f"GFF3      : {args.gff}")
    logging.info(f"Flank     : {args.flank} bp")
    logging.info(f"Gene scope: {'all genes' if args.all_genes else 'sorgoleone targets'}\n")

    dmr_df = pd.read_csv(args.dmr, sep="\t")
    logging.info(f"Loaded {len(dmr_df)} DMRs from {args.dmr}\n")

    loc_names   = None if args.all_genes else LOC_NAMES
    label_to_egi = None if args.all_genes else LABEL_TO_EGI

    feature_dir = outdir / "features"
    feature_beds = build_feature_beds(args.gff, loc_names, args.flank, feature_dir)

    logging.info(f"\nIntersecting {len(dmr_df)} DMRs against {len(feature_beds)} feature types...")
    dmr_annotated = annotate_dmrs(dmr_df, feature_beds, outdir, label_to_egi=label_to_egi)

    print_summary(dmr_annotated)

    out_path = Path(args.out) if args.out else outdir / "DMR_annotated.tsv"
    dmr_annotated.to_csv(out_path, sep="\t", index=False)
    logging.info(f"\nAnnotated DMR table saved to: {out_path}")


if __name__ == "__main__":
    main()
