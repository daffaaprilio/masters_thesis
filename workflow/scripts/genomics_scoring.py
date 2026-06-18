#!/opt/conda/bin/python3

import numpy as np
import cyvcf2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Reads VCF file, assigns scores to each gene containing variants listed in that VCF file.
        Variants are associated with the samples' phenotype, i.e., SBC10-private variants are associated with low TAA production (aconitate isomerase), SBC11-private variants; D gene (stem juiciness).
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input",
        help="Path to the input sample-exclusive VCF file",
    )

    parser.add_argument(
        "-o", "--output", 
        help="Parent directory of the output files",
    )

# prepare dataframe
ANN_FIELDS = [
    "ann_allele", "ann_effect", "ann_impact", "ann_gene_name", "ann_gene_id",
    "ann_feature_type", "ann_feature_id", "ann_biotype", "ann_rank",
    "ann_hgvs_c", "ann_hgvs_p", "ann_cdna_pos", "ann_cds_pos", "ann_aa_pos",
    "ann_distance", "ann_extra",
]

SIFT_FIELDS = [
    "sift_allele", "sift_transcript", "sift_gene_id", "sift_gene_name",
    "sift_region", "sift_variant_type", "sift_aa_change", "sift_aa_pos",
    "sift_score", "sift_median", "sift_num_seqs", "sift_allele_type",
    "sift_prediction",
]

def _parse_ann(raw: str) -> dict:
    records = []
    for entry in raw.split(","): # a row could have one or more ANN annotation, because ANN annotates every transcripts affected by that variant site
        parts = entry.split("|")
        parts += [""] * (len(ANN_FIELDS) - len(parts))
        records.append(dict(zip(ANN_FIELDS, parts[:len(ANN_FIELDS)])))
    return records

def _parse_siftinfo(raw:str) -> dict:
    