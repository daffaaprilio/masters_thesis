#!/Users/daffa/miniconda3/envs/sbi/bin/python3

import cyvcf2
import pandas as pd
import argparse
import logging
from datetime import datetime
from pathlib import Path

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]

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

    stem = Path(args.vcf).name.split(".vcf")[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = Path(__file__).parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"parse_vcf_to_tsv.{stem}.{timestamp}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    log = logging.getLogger(__name__)
    log.info("Input:  %s", args.vcf)
    log.info("Outdir: %s", args.outdir)

    vcf = cyvcf2.VCF(args.vcf, samples=SAMPLES)

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
            # Genotypes: list of [allele1, allele2, phased] per sample
            "GT_SBC4": variant.genotypes[0],
            "GT_SBC10": variant.genotypes[1],
            "GT_SBC11":  variant.genotypes[2],
            "GT_SBC23": variant.genotypes[3],
        }

        for ann in annotations:
            rows.append({**base, **ann})

    out_path = Path(args.outdir) / f"{stem}.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    log.info("Saved %d rows to %s", len(rows), out_path)
    log.info("Log written to %s", log_file)

if __name__ == "__main__":
    main()