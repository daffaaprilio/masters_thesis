#!/Users/daffa/miniconda3/envs/sbi/bin/python3

# Convert SnpEff gene IDs (LocusTag e.g. SORBI_3001G000100) in the ANN INFO field
# to NCBI GeneIDs using the filtered gene_info TSV as a lookup table.

import argparse
import csv
import sys


def load_locus_to_geneid(gene_info_path: str) -> dict:
    """Return {LocusTag: GeneID} from a filtered NCBI gene_info TSV."""
    mapping = {}
    with open(gene_info_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            locus = row.get("LocusTag", "").strip()
            gene_id = row.get("GeneID", "").strip()
            if locus and locus != "-":
                mapping[locus] = gene_id
    return mapping


def convert_ann_field(ann_value: str, mapping: dict) -> str:
    """Replace gene name (index 3) and gene ID (index 4) in every ANN sub-entry."""
    converted = []
    for entry in ann_value.split(","):
        parts = entry.split("|")
        if len(parts) > 4:
            locus = parts[3]
            gene_id = mapping.get(locus, locus)  # fall back to original if not found
            parts[3] = gene_id
            parts[4] = gene_id
        converted.append("|".join(parts))
    return ",".join(converted)


def process_info(info: str, mapping: dict) -> str:
    """Rewrite only the ANN= sub-field inside the INFO column."""
    fields = info.split(";")
    result = []
    for field in fields:
        if field.startswith("ANN="):
            field = "ANN=" + convert_ann_field(field[4:], mapping)
        result.append(field)
    return ";".join(result)


def main():
    parser = argparse.ArgumentParser(
        description="Convert SnpEff LocusTags in ANN field to NCBI GeneIDs."
    )
    parser.add_argument("--input",  required=True, help="SnpEff-annotated VCF file")
    parser.add_argument("--output", required=True, help="Output VCF file path")
    parser.add_argument("--key",    required=True, help="Filtered NCBI gene_info TSV (output of filter_gene_info rule)")
    args = parser.parse_args()

    mapping = load_locus_to_geneid(args.key)
    print(f"Loaded {len(mapping):,} LocusTag → GeneID entries", file=sys.stderr)

    converted = 0
    with open(args.input) as vcf_in, open(args.output, "w") as vcf_out:
        for line in vcf_in:
            if line.startswith("#"):
                vcf_out.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 7:
                cols[7] = process_info(cols[7], mapping)
                converted += 1
            vcf_out.write("\t".join(cols) + "\n")

    print(f"Done — converted {converted:,} variant lines → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
