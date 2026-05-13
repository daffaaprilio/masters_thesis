# Exploration of Sorgoleone Key Genes
Key gene sequence information: [Supplementary information](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Ftpj.16263&file=tpj16263-sup-0002-TableS1-S3.xlsx) of [Maharjan et al. 2023](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16263). <br>
Sorgoleone biosynthesis pathway key genes/enzymes, obtained from BLAST-ing gene sequences:
1.  SbDES2 <br>
    [Associated gene](https://www.ncbi.nlm.nih.gov/gene?term=EF206347[Nucleotide%20Accession]&RID=05AT5F01016&log$=genealign&blast_rank=1): LOC8066368 aka SORBI_3004G260600
2.  SbDES3 <br>
    [Associated gene](https://www.ncbi.nlm.nih.gov/gene?term=NM_001424099[Nucleotide%20Accession]&RID=05B1Z21G014&log$=genealign&blast_rank=1): LOC8079957 aka SORBI_3005G002700
3.  SbARS1 <br>
    No associated gene found
4.  SbOMT3 <br>
    [Associated gene](https://www.ncbi.nlm.nih.gov/gene?term=NM_001424100[Nucleotide%20Accession]&RID=05BEKV90014&log$=genealign&blast_rank=1): LOC8080259 aka SORBI_3006G007900
5.  SbCYP71AM1 <br>
    [Associated gene](https://www.ncbi.nlm.nih.gov/gene?term=XM_002451987[Nucleotide%20Accession]&RID=05BH80V7014&log$=genealign&blast_rank=2): LOC8081692 aka SORBI_3004G139300

## Comparative genomics analysis
Technical steps
1.  `03_merge_vcf.sh`: Merging VCF files for all samples, then filter variant sites to only include the following genes. <br>
    ```shell
    ./analysis/03_merge_vcf.sh # merges phased VCF from all 4 samples
    ```
2.  `01_annotate_vcf.sh`: Annotate (and rename) merged VCF (chromosome names have to be renamed beforehand, so that SnpEFF can correctly annotate them. Yes, SnpEFF database has different chromosome naming system.) <br>
    ```shell
    ./analysis/01_annotate_vcf.sh merged /Users/daffa/workspace/infobio/thesis/analysis/data/vcf/merged/merged.phased.vcf.gz
    ```
3.  `04_annot_vcf_to_tsv.py`: Parse VCF into TSV to explore in a notebook environment
    ```shell
    python3 /Users/daffa/workspace/infobio/thesis/analysis/04_annot_vcf_to_tsv.py \
    -v /Users/daffa/workspace/infobio/thesis/analysis/data/vcf/annotated/merged.annotated.vcf.gz \
    -o /Users/daffa/workspace/infobio/thesis/analysis/data/tsv
    ```
4.  `sorgoleone/sorgoleone.ipynb`: Analyze in notebook: first, filter variant sites from sorgholeone genes only. <br>
    Group variant sites by [impact, i.e., MODERATE, LOW, MODIFIER, etc.](https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files)

Biological insights
-   **Frameshift mutation in SbDES2 gene found exclusively in SBC4.** <br>
    Location (bp): 4:60,583,662 (in IGV NC_012873.2:60,583,642-60,583,681)
