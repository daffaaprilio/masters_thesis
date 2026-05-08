#!/usr/bin/python3

import argparse
import glob
import os
import shutil
import tempfile
from pyspark.sql import SparkSession
from pyspark.sql.functions import col

# NCBI All_Data.gene_info column order (header line begins with '#tax_id')
GENE_INFO_COLUMNS = [
    "tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs",
    "chromosome", "map_location", "description", "type_of_gene",
    "Symbol_from_nomenclature_authority", "Full_name_from_nomenclature_authority",
    "Nomenclature_status", "Other_designations", "Modification_date", "Feature_type",
]

SORGHUM_BICOLOR_TAXID = 4558


def build_spark(driver_memory: str = "12g") -> SparkSession:
    return (
        SparkSession.builder
        .appName("extract_sbi_ncbi_gene_info")
        .master("local[*]")
        .config("spark.driver.memory", driver_memory)
        # avoid writing _SUCCESS / .crc sidecar files
        .config("spark.hadoop.mapreduce.fileoutputcommitter.marksuccessfuljobs", "false")
        .getOrCreate()
    )


def main():
    parser = argparse.ArgumentParser(
        description="Filter NCBI gene_info for Sorghum bicolor (taxon 4558)."
    )
    parser.add_argument("--input", required=True, help="Path to All_Data.gene_info (TSV, ~9 GB)")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--tax-id", type=int, default=SORGHUM_BICOLOR_TAXID,
        help=f"NCBI taxon ID to keep (default: {SORGHUM_BICOLOR_TAXID})",
    )
    parser.add_argument(
        "--driver-memory", default="12g",
        help="Spark driver memory (default: 12g)",
    )
    args = parser.parse_args()

    spark = build_spark(args.driver_memory)
    spark.sparkContext.setLogLevel("WARN")

    # gene_info header starts with '#tax_id' — read as plain text, skip comment,
    # then split on tab to get clean column names.
    raw = spark.read.text(args.input)

    header_row = raw.filter(col("value").startswith("#tax_id")).limit(1)
    data_rows  = raw.filter(~col("value").startswith("#"))

    # Parse header to derive column positions at runtime (guards against future
    # NCBI schema additions).
    header_str = header_row.collect()[0]["value"]
    columns = [c.lstrip("#") for c in header_str.split("\t")]

    from pyspark.sql.functions import split as spark_split
    split_col = spark_split(col("value"), "\t")
    df = data_rows.select(
        *[split_col.getItem(i).alias(columns[i]) for i in range(len(columns))]
    )

    filtered = df.filter(col("tax_id") == str(args.tax_id))

    tmp_dir = tempfile.mkdtemp()
    try:
        (
            filtered.coalesce(1).write
            .option("header", "true")
            .option("sep", "\t")
            .option("emptyValue", "-")
            .mode("overwrite")
            .csv(tmp_dir)
        )
        part_file = glob.glob(os.path.join(tmp_dir, "part-*.csv"))[0]
        shutil.move(part_file, args.output)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    count = filtered.count()
    print(f"Done — wrote {count:,} rows for taxon {args.tax_id} to {args.output}")

    spark.stop()


if __name__ == "__main__":
    main()
