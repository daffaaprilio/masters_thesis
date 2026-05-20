# Identification of TAA-related genes
Functional genomics and epigenomics, both enhanced with gene co-expression analysis to identify key genes in the biosynthesis pathway as well as secretion of trans-aconitic acid (TAA) in sorghum.

## Comparative genomics analysis
Technical steps
1.  `snpeff_prep.sh`: Run the script to prepare the SnpEff database first <br>
    ```shell
    ./docker/run.sh bash analysis/scripts/snpeff_prep.sh
    ```

2.  `private_variants.sh` Find variant sites private to each sample. for each file in the processed VCF files directory (`results/vcf_processing/SBC*.phased.vcf.gz`). <br>
    ```shell
    ./docker/run.sh bash analysis/scripts/private_variants.sh
    ```

3.  Annotate private variants using SnpEFF. <br>
    ```shell
    ./docker/run.sh analysis/scripts/annotate_vcf.sh SBC10.private analysis/data/vcf/private_variants/SBC10.private.vcf.gz analysis/data/vcf/private_variants
    ```

4.  `annot_single_vcf_to_tsv.py`: Parse single-sampled VCF into TSV to explore in a notebook environment. The other file is designed to parse multi-sample (4 samples merged) VCF. <br>
    ```shell
    ./docker/run.sh python3 analysis/scripts/annot_single_vcf_to_tsv.py -v analysis/data/vcf/private_variants/SBC10.private.annotated.vcf.gz -o analysis/data/tsv
    ```

5.  `03_TAA/TAA.ipynb`: Analyze in notebook.