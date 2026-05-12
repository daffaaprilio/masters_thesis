#!/usr/bin/env python3
import argparse

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

ANN_FIELDS = [
    "Allele", "Effect", "Impact", "Gene_Name", "Gene_ID",
    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
    "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "AA_pos",
    "Distance", "Errors",
]

OUT_COLS = [
    "chrom", "pos", "ref", "alt", "qual",
    "gene", "effect", "impact", "biotype", "hgvs_c", "hgvs_p",
    "GT", "GQ", "DP", "AF",
]


def best_ann(ann_raw):
    """Return the ANN entry with highest impact that is protein_coding HIGH/MODERATE."""
    candidates = []
    for entry in ann_raw.split(","):
        parts = entry.split("|")
        parts += [""] * (len(ANN_FIELDS) - len(parts))
        ann = dict(zip(ANN_FIELDS, parts[: len(ANN_FIELDS)]))
        if ann["Transcript_BioType"] == "protein_coding" and ann["Impact"] in ("HIGH", "MODERATE"):
            candidates.append(ann)
    if not candidates:
        return None
    candidates.sort(key=lambda a: IMPACT_RANK.get(a["Impact"], 99))
    return candidates[0]


def parse_format(fmt_str, sample_str):
    return dict(zip(fmt_str.split(":"), sample_str.split(":")))


def main():
    ap = argparse.ArgumentParser(description="Filter annotated VCF for high-confidence candidate variants.")
    ap.add_argument("--input",  required=True, help="Annotated EGI VCF (uncompressed)")
    ap.add_argument("--output", required=True, help="Output TSV")
    ap.add_argument("--sample", default=None,  help="Sample column name (default: first sample)")
    args = ap.parse_args()

    vcf_cols = None
    sample_col = args.sample
    records = []

    with open(args.input) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                vcf_cols = line.lstrip("#").split("\t")
                if sample_col is None:
                    sample_col = vcf_cols[9]
                continue

            fields = dict(zip(vcf_cols, line.split("\t")))

            # 1. QUAL >= 20
            try:
                qual = float(fields["QUAL"])
            except (ValueError, KeyError):
                continue
            if qual < 20:
                continue

            fmt = parse_format(fields["FORMAT"], fields[sample_col])

            # 2. 10 <= DP <= 100
            try:
                dp = int(fmt["DP"])
            except (ValueError, KeyError):
                continue
            if not (10 <= dp <= 100):
                continue

            # 5. GT == 1/1 (accept phased 1|1 too)
            gt = fmt.get("GT", "")
            if gt.replace("|", "/") != "1/1":
                continue

            # 6. AF >= 0.75
            try:
                af = float(fmt["AF"])
            except (ValueError, KeyError):
                continue
            if af < 0.75:
                continue

            # 3 & 4. protein_coding + HIGH/MODERATE impact (via best_ann)
            ann_raw = None
            for tok in fields["INFO"].split(";"):
                if tok.startswith("ANN="):
                    ann_raw = tok[4:]
                    break
            if ann_raw is None:
                continue

            ann = best_ann(ann_raw)
            if ann is None:
                continue

            records.append({
                "chrom":   fields["CHROM"],
                "pos":     int(fields["POS"]),
                "ref":     fields["REF"],
                "alt":     fields["ALT"],
                "qual":    qual,
                "gene":    ann["Gene_Name"],
                "effect":  ann["Effect"],
                "impact":  ann["Impact"],
                "biotype": ann["Transcript_BioType"],
                "hgvs_c":  ann["HGVS_c"],
                "hgvs_p":  ann["HGVS_p"],
                "GT":      fmt.get("GT", ""),
                "GQ":      fmt.get("GQ", ""),
                "DP":      dp,
                "AF":      af,
            })

    records.sort(key=lambda r: (IMPACT_RANK.get(r["impact"], 99), r["gene"], r["pos"]))

    print(f"Passing variants: {len(records)}", flush=True)

    with open(args.output, "w") as out:
        out.write("\t".join(OUT_COLS) + "\n")
        for r in records:
            out.write("\t".join(str(r[c]) for c in OUT_COLS) + "\n")


if __name__ == "__main__":
    main()
