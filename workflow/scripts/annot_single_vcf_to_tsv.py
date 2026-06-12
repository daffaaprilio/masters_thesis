#!/usr/bin/env python3

import cyvcf2
import pandas as pd
import argparse
import logging
from pathlib import Path

ANN_FIELDS = [
    "allele", "effect", "impact", "gene_name", "gene_id",
    "feature_type", "feature_id", "biotype", "rank",
    "hgvs_c", "hgvs_p", "cdna_pos", "cds_pos", "aa_pos",
    "distance", "extra"
] # from https://pcingola.github.io/SnpEff/snpeff/inputoutput/#ann-field-vcf-output-files

def parse_ann(ann_string):
    """Parse one ANN= value into a list of dicts (one per annotation)."""
    records = []
    for annotation in ann_string.split(","):
        parts = annotation.split("|")
        # Pad to 16 fields if trailing pipes are missing
        parts += [""] * (16 - len(parts))
        records.append(dict(zip(ANN_FIELDS, parts[:16])))
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcf", help="Annotated VCF file to convert")
    parser.add_argument("-o", "--outdir", help="Path of the output file")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger(__name__)
    log.info("Input:  %s", args.vcf)
    log.info("Outdir: %s", args.outdir)

    vcf = cyvcf2.VCF(args.vcf)
    samples = vcf.samples
    log.info("Samples in VCF: %s", samples)

    rows = []
    for variant in vcf:
        ann_raw = variant.INFO.get("ANN")
        annotations = parse_ann(ann_raw) if ann_raw else [{}]

        base = {
            "chrom": variant.CHROM,
            "pos":   variant.POS,
            "ref":   variant.REF,
            "alt":   ",".join(variant.ALT),
            "qual":  variant.QUAL,
        }
        for i, sample in enumerate(samples):
            a1, a2, phased = variant.genotypes[i]
            sep = "|" if phased else "/"
            a1_str = "." if a1 < 0 else str(a1)
            a2_str = "." if a2 < 0 else str(a2)
            base[f"GT_{sample}"] = f"{a1_str}{sep}{a2_str}"

        for ann in annotations:
            rows.append({**base, **ann})

    out_path = Path(args.outdir) / f"{stem}.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    log.info("Saved %d rows to %s", len(rows), out_path)

if __name__ == "__main__":
    main()