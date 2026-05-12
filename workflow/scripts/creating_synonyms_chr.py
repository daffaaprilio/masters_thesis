#!/Users/daffa/miniconda3/envs/sbi/bin/python3

import pandas as pd

WDIR = "/Users/daffa/workspace/infobio/thesis"
vcf_list_path = f"{WDIR}/workflow/scripts/vcf_chr_list.txt"
snpeff_db_list_path = f"{WDIR}/workflow/scripts/snpeff_db_chr_list.txt"
output_synonym_path = f"{WDIR}/workflow/scripts/synonyms.txt"

def list_vcf_chr(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("##contig="):
                parts = line.strip().rstrip(">").split(",")
                contig_id = parts[0].replace("##contig=<ID=", "")
                length = int(parts[1].replace("length=", ""))
                rows.append({"contig": contig_id, "length": length})
    return pd.DataFrame(rows, columns=["contig", "length"])

def list_snpeff_db_chr(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4 and parts[0] == "#" and parts[3] == "Standard":
                contig_id = parts[1].strip("'")
                length = int(parts[2])
                rows.append({"contig": contig_id, "length": length})
    return pd.DataFrame(rows, columns=["contig", "length"])

df_vcf = list_vcf_chr(vcf_list_path)
df_snpeff = list_snpeff_db_chr(snpeff_db_list_path)

df_synonyms = pd.merge(
    df_vcf, df_snpeff,
    on="length",
    suffixes=("_vcf", "_snpeff")
)

df_synonyms[["contig_vcf", "contig_snpeff"]].to_csv(output_synonym_path, index=False, header=False, sep="\t")

