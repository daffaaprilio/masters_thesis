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

## Key genes and their homologs

| Key Gene | Gene (EGI) | Locus name (SORBI_) | Chr : Location (bp) | Strand |
|---|---|---|---|---|
| **SbDES2** | 8066368 | SORBI_3004G260600 | Chr4 NC_012873.2: 60,583,460–60,584,984 | + |
| ↳ *SbDES2-homolog* | 110435045 | SORBI_3004G260800 | Chr4 NC_012873.2: 60,620,793–60,625,730 | + |
| **SbDES3** | 8079957 | SORBI_3005G002700 | Chr5 NC_012874.2: 207,074–209,472 | − |
| ↳ *SbDES3-homolog-1* | 8072903 | SORBI_3008G002800 | Chr8 NC_012877.2: 248,854–251,593 | − |
| ↳ *SbDES3-homolog-2* | 8055482 | SORBI_3008G003200 | Chr8 NC_012877.2: 280,238–285,184 | − |
| ↳ *SbDES3-homolog-3* | 8079958 | SORBI_3005G002800 | Chr5 NC_012874.2: 218,355–222,920 | − |
| **SbOMT3** | 8080259 | SORBI_3006G007900 | Chr6 NC_012875.2: 1,183,885–1,185,775 | − |
| ↳ *SbOMT3-homolog-1* | 8076922 | SORBI_3005G086600 | Chr5 NC_012874.2: 12,031,041–12,032,827 | − |
| ↳ *SbOMT3-homolog-2* | 110436225 | SORBI_3006G008000 | Chr6 NC_012875.2: 1,203,224–1,205,498 | − |
| ↳ *SbOMT3-homolog-3* | 8085153 | SORBI_3007G074800 | Chr7 NC_012876.2: 8,583,133–8,585,000 | + |
| **SbCYP71AM1** | 8081692 | SORBI_3004G139300 | Chr4 NC_012873.2: 39,976,273–39,978,067 | + |

## Comparative genomics analysis
Technical steps
1.  `snpeff_prep.sh`: Run the script to prepare the SnpEff database first <br>
    ```shell
    ./docker/run.sh bash analysis/scripts/snpeff_prep.sh
    ```

2.  `merge_vcf.sh`: Merging VCF files for all samples, then filter variant sites to only include the following genes. <br>
    ```shell
    ./docker/run.sh bash analysis/scripts/merge_vcf.sh
    ```

3.  `annotate_vcf.sh`: Annotate (and rename) merged VCF (chromosome names have to be renamed beforehand, so that SnpEFF can correctly annotate them. Yes, SnpEFF database has different chromosome naming system.) <br>
    ```shell
    ./docker/run.sh bash analysis/scripts/annotate_vcf.sh \
        merged \
        analysis/data/vcf/merged/merged.phased.vcf.gz \
        analysis/data/vcf/annotated
    ```

4.  `annot_vcf_to_tsv.py`: Parse VCF into TSV to explore in a notebook environment
    ```shell
    ./docker/run.sh python3 analysis/scripts/annot_vcf_to_tsv.py \
        -v analysis/data/vcf/annotated/merged.annotated.vcf.gz \
        -o analysis/data/tsv
    ```

5.  `04_sorgoleone/sorgoleone.ipynb`: Analyze in notebook: first, filter variant sites from sorgholeone genes only. <br>
    Group variant sites by [impact, i.e., MODERATE, LOW, MODIFIER, etc.](https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files)

## Comparative epigenomics analysis
Technical steps

1.  `filter_bedmethyl_sorgoleone.py`: Filter bedMethyl files to regions corresponding to the 4 sorgoleone pathway genes (SbDES2, SbDES3, SbOMT3, SbCYP71AM1) and their 7 co-expressed homologs, each extended by a 2000 bp upstream flank (strand-aware) to capture regulatory regions, i.e., promoter. Gene coordinates are extracted from the GFF3 annotation, then `bedtools intersect` is applied to each sample's filtered bedMethyl file. <br>
    ```shell
    ./docker/run.sh python analysis/scripts/filter_bedmethyl_sorgoleone.py \
        --gff resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff \
        --bedmethyl resources/bedmethyl/SBC4.filtered.bed \
                    resources/bedmethyl/SBC10.filtered.bed \
                    resources/bedmethyl/SBC11.filtered.bed \
                    resources/bedmethyl/SBC23.filtered.bed \
        --outdir analysis/data/sorgoleone_bedmethyl
    ```

2.  `validate_bedmethyl_sorgoleone.py`: Another check on the methylation status. Check saved in `analysis/data/sorgoleone_bedmethyl/sorgoleone_bedmethyl_validation.tsv` <br>
    ```shell
    ./docker/run.sh python analysis/scripts/validate_bedmethyl_sorgoleone.py \
        --gff resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff
    ```

3.  `bedmethyl_to_dss.py`: Convert filtered bedMethyl files to DSS input format. For 5mC, collapses Watson/Crick CpG strand pairs by summing coverage and modified read counts. For 6mA, keeps strands separate. Excludes artefact code `21839` and sites with `valid_coverage < 5`. Outputs `{sample}.5mC.dss.txt` and `{sample}.6mA.dss.txt`. <br>
    ```shell
    ./docker/run.sh python analysis/scripts/bedmethyl_to_dss.py
    ```

4.  `run_dss_dmr.R`: Run pairwise DSS DMR analysis across all six accession pairs (SBC04 vs SBC10, SBC04 vs SBC11, SBC04 vs SBC23, SBC10 vs SBC11, SBC10 vs SBC23, SBC11 vs SBC23). Each pair is tested with smoothed DML testing (`smoothing.span=500`) followed by DMR calling (`delta=0.1`, `p<0.2`, `minlen=50`, `minCG=20`, `dis.merge=300`). Analysis is exploratory (n=1 per group, no replicates). Outputs one DMR and one DML TSV per pair, plus a combined DMR table and summary. <br>
    ```shell
    ./docker/run.sh Rscript analysis/scripts/run_dss_dmr.R
    ```

5.  `annotate_DMR.py`: Annotate each DMR with the genomic feature(s) it overlaps — promoter (2 kb strand-aware upstream flank), exon, CDS, or intron. Feature BEDs are derived from the GFF3 via parent-child traversal (`gene → mRNA → exon/CDS`) scoped to the 11 target genes. Introns are computed per gene as gene body minus merged exons (`bedtools subtract`). UTR features are absent from this NCBI annotation. Outputs `DMR_annotated.tsv` with two added columns: `features` (semicolon-separated overlap types, or `intergenic`) and `gene_label`. <br>
    ```shell
    ./docker/run.sh python analysis/scripts/annotate_DMR.py
    ```


## Gene co-expression analysis
1.  Listing genes (ready to copy to ATTED-II coexpression network illustrator)
    ```
    8066368
    8079957
    8080259
    8081692

    # extra (homologs genes), but the visualized networks cluster similarly, though
    110435045
    8072903
    8055482
    8079958
    8076922
    110436225
    8085153
    ```
2.  Draw co-expression network using [ATTED-II NetworkDrawer](https://atted.jp/top_draw/#NetworkDrawer)

## Biological insights

-   **Gene coexpression calculation returns clusters of highly co-expressed homologous loci**
    - We hypothesized that the four identified key sorgoleone biosnynthesis pathway genes are highly co-expressed (add more proofs, partially its physical interaction in the cell (Maharjan et al. (2023) Interaction between sorgoleone biosynthesis enzymes section, last paragraph)).
    - However, it seems to be no direct co-expression between those 4 genes. Instead, 3 out of 4 genes show high co-expression z-score within their homologs (shown in the next two sections).

-   **Functional disruption of SbDES2 and predicted compensatory mechanism in SBC4.** <br>
    - Comparative genomic analysis revealed a critical structural variant in the SbDES2 gene (SORBI_3004G260600, gene ID: 8066368) within the SBC4 accession. Specifically, a single-nucleotide insertion of guanine at position 60,583,662 on chromosome 4 (Location (chr:bp): NC_012873.2:60,583,662) resulted in a frameshift mutation that is predicted to disrupt normal protein function (SnpEFF prediction). 
        - Further prediction, particularly to answer the question whether the mutation actually disrupts the gene's sorgoleone biosynthesis pathway function, is necessary. A good idea might be to check the predicted protein function (i.e., if the mutation happens at the critical domain that catalyzes the fatty acid desaturation function).
    - Gene co-expression network analysis identified a functionally homologous gene (gene ID: 110435045, blastp identity 89.38) that exhibits similar expression patterns to SbDES2 (coexpression z-score: 16.4). Importantly, variant analysis of this homologous locus revealed no high-impact mutations (SnpEFF annotation), suggesting that the protein product remains structurally and functionally intact in SBC4.
    - Further protein domain prediction, followed by an experimental validation will be necessary to confirm the predicted loss-of-function in SbDES2 (SORBI_3004G260600, gene ID: 8066368) and demonstrate functional compensation by its homolog (SORBI_3004G260800, gene ID: 110435045).

-   **Lost of start codon in the SbDES3 homologous locus ("8079957": "SORBI_3005G002700")** <br>
    - Gene co-expression calculation reveals that SbDES3 is highly co-expressed with the other 3 genes, which are homologous to each other.
        - 8072903 aka SORBI_3008G002800, blastp identity 90.81, coex z-score: 15.5, probable FA desaturase DES1
        - 8055482 aka SORBI_3008G003200, blastp identity 84.22, coex z-score: 15.5, FA desaturase DES3-like
        - 8079958 aka SORBI_3005G002800, blastp identity 86.56, coex z-score: 15.5, FA desaturase DES3-like
    - A critical SNP was found in 8055482 in the SBC4 sample (chromosome 8, NC_012877.2:284,923, which is a start codon negative direction). A mutated into G, causing the lost of start codon.

-   **Insertion and deletion in SbOMT3-like loci**
    - Gene co-expression network suggests that SbOMT3 (SORBI_3006G007900, gene ID: 8080259) is highly co-expressed with 3 other genes, which are homologous to each other
        - 8076922 aka SORBI_3005G086600, blastp identity 92.51, coex z-score: 8.5, 5-pentadecatrienyl resorcinol O-methyltransferase-like
        - 110436225 aka SORBI_3006G008000, blastp identity 97.86, coex z-score: 8.5, 5-pentadecatrienyl resorcinol O-methyltransferase-like
        - 8085153 aka SORBI_3007G074800, blastp identity 70, coex z-score: 8.5, 5-pentadecatrienyl resorcinol O-methyltransferase-like
    - A 50 bp-long deletion was found in 110436225 (a heterozygous variant) in SBC11 (NC_012875.2:1,205,067-1,205,117), causing a frameshift mutation
    - Insertion of a C base in NC_012876.2:8,584,480 was observed in SBC4 and SBC10. This insertion happens in the splice acceptor region and annotated as highly impactful.

-   **Hypermethylated regions in the promoter of SbDES2-homolog (gene ID 110435045) of SBC10 and SBC11**
    - SBC10 & SBC11's SbDES2-homolog gene's promoter region is shown to be hypermethylated, compared to SBC23. This suggests that either the expression of the SbDES2-homolog gene (gene ID 110435045) in SBC10 & SBC11 is suppressed or the expression of that gene is high in SBC23.


-   **Hypermethylated region in the promoter of SbOMT3-homolog (gene ID 8085153) of SBC10**
    - Hypermethylation in the promoter region of the SbOMT3-homolog gene (gene ID 8085153) in SBC10, compared to SBC23, suggests that either the expression of one (gene ID 8085153) of the SbOMT3-homolog genes in SBC10 is suppressed or the expression of that gene is high in SBC23.

### Summary

**SBC4**
| Gene | Gene ID | Locus Name | High-impact variants | DMR |
|------|---------|------------|----------------------|-----|
| SbDES2 | 8066368 | SORBI_3004G260600 | Frameshift — insertion G at NC_012873.2:60,583,662 | — |
| SbDES2-homolog | 110435045 | SORBI_3004G260800 | — | — |
| SbDES3 | 8079957 | SORBI_3005G002700 | — | — |
| SbDES3-homolog-1 | 8072903 | SORBI_3008G002800 | — | — |
| SbDES3-homolog-2 | 8055482 | SORBI_3008G003200 | Start codon loss — A→G at NC_012877.2:284,923 | — |
| SbDES3-homolog-3 | 8079958 | SORBI_3005G002800 | — | — |
| SbOMT3 | 8080259 | SORBI_3006G007900 | — | — |
| SbOMT3-homolog-1 | 8076922 | SORBI_3005G086600 | — | — |
| SbOMT3-homolog-2 | 110436225 | SORBI_3006G008000 | — | — |
| SbOMT3-homolog-3 | 8085153 | SORBI_3007G074800 | Splice acceptor insertion — C at NC_012876.2:8,584,480 | — |
| SbCYP71AM1 | 8081692 | SORBI_3004G139300 | — | — |

**SBC10**
| Gene | Gene ID | Locus Name | High-impact variants | DMR |
|------|---------|------------|----------------------|-----|
| SbDES2 | 8066368 | SORBI_3004G260600 | — | — |
| SbDES2-homolog | 110435045 | SORBI_3004G260800 | — | Hypermethylated promoter vs SBC23 |
| SbDES3 | 8079957 | SORBI_3005G002700 | — | — |
| SbDES3-homolog-1 | 8072903 | SORBI_3008G002800 | — | — |
| SbDES3-homolog-2 | 8055482 | SORBI_3008G003200 | — | — |
| SbDES3-homolog-3 | 8079958 | SORBI_3005G002800 | — | — |
| SbOMT3 | 8080259 | SORBI_3006G007900 | — | — |
| SbOMT3-homolog-1 | 8076922 | SORBI_3005G086600 | — | — |
| SbOMT3-homolog-2 | 110436225 | SORBI_3006G008000 | — | — |
| SbOMT3-homolog-3 | 8085153 | SORBI_3007G074800 | Splice acceptor insertion — C at NC_012876.2:8,584,480 | Hypermethylated promoter vs SBC23 |
| SbCYP71AM1 | 8081692 | SORBI_3004G139300 | — | — |

**SBC11**
| Gene | Gene ID | Locus Name | High-impact variants | DMR |
|------|---------|------------|----------------------|-----|
| SbDES2 | 8066368 | SORBI_3004G260600 | — | — |
| SbDES2-homolog | 110435045 | SORBI_3004G260800 | — | Hypermethylated promoter vs SBC23 |
| SbDES3 | 8079957 | SORBI_3005G002700 | — | — |
| SbDES3-homolog-1 | 8072903 | SORBI_3008G002800 | — | — |
| SbDES3-homolog-2 | 8055482 | SORBI_3008G003200 | — | — |
| SbDES3-homolog-3 | 8079958 | SORBI_3005G002800 | — | — |
| SbOMT3 | 8080259 | SORBI_3006G007900 | — | — |
| SbOMT3-homolog-1 | 8076922 | SORBI_3005G086600 | — | — |
| SbOMT3-homolog-2 | 110436225 | SORBI_3006G008000 | Frameshift — 50 bp deletion at NC_012875.2:1,205,067–1,205,117 | — |
| SbOMT3-homolog-3 | 8085153 | SORBI_3007G074800 | — | — |
| SbCYP71AM1 | 8081692 | SORBI_3004G139300 | — | — |

**SBC23**
| Gene | Gene ID | Locus Name | High-impact variants | DMR |
|------|---------|------------|----------------------|-----|
| SbDES2 | 8066368 | SORBI_3004G260600 | — | — |
| SbDES2-homolog | 110435045 | SORBI_3004G260800 | — | — |
| SbDES3 | 8079957 | SORBI_3005G002700 | — | — |
| SbDES3-homolog-1 | 8072903 | SORBI_3008G002800 | — | — |
| SbDES3-homolog-2 | 8055482 | SORBI_3008G003200 | — | — |
| SbDES3-homolog-3 | 8079958 | SORBI_3005G002800 | — | — |
| SbOMT3 | 8080259 | SORBI_3006G007900 | — | — |
| SbOMT3-homolog-1 | 8076922 | SORBI_3005G086600 | — | — |
| SbOMT3-homolog-2 | 110436225 | SORBI_3006G008000 | — | — |
| SbOMT3-homolog-3 | 8085153 | SORBI_3007G074800 | — | — |
| SbCYP71AM1 | 8081692 | SORBI_3004G139300 | — | — |