# Identification of TAA-related genes
Functional genomics and epigenomics, both enhanced with gene co-expression analysis to identify key genes in the biosynthesis pathway as well as secretion of trans-aconitic acid (TAA) in sorghum.

## Sample Phenotypes

| Sample | TAA Production | TAA Secretion | Callus Formation |
|--------|---------------|---------------|-----------------|
| SBC4 | ++ | High | Mid |
| SBC10 | +++ | Low | Good |
| SBC11 | - | High | Mid |
| SBC23 | ++ | High | Good |

## Study Design
### Data Preparation
```mermaid
flowchart TD
    subgraph ACC["Four sorghum accessions"]
        SBC04["SBC04"]
        SBC10["SBC10"]
        SBC11["SBC11"]
        SBC23["SBC23"]
    end

    SBC10 -->|"ONT R10.4.1 duplex"| SEQ["Long-read sequencing\n(Dorado v1.3.0, SUP model)"]
    SBC04 -->|"ONT R10.4.1 duplex"| SEQ
    SBC11 -->|"ONT R10.4.1 duplex"| SEQ
    SBC23 -->|"ONT R10.4.1 duplex"| SEQ

    SEQ -->|"BTx623 reference genome"| VAR["Variant calling\n(Clair3)"]
    SEQ --> METH["Methylation profiling\n(modkit pileup)"]

    METH --> DMR["Differential methylation analysis\n(R DSS package)"]
    VAR --> VCFP["VCF Processing\n(i.e., filtering, WhatsHap phasing)"]
    VCFP --> ISEC["Sample-private variant profiling\n(bcftools isec)"]

    subgraph SBC10RES["SBC10 analysis outputs"]
        SBC10DMR["Differential methylation\nSBC10 vs others"]
        VARPSBC10["SBC10-private variants"]
    end

    DMR --> SBC10DMR
    ISEC --> VARPSBC10

    style ACC fill:#e6e6e6
```

### Analysis Overview
```mermaid
flowchart TD
 subgraph s1["Genomics approach"]
        VARPSBC10["SBC10-private variants"]
        ANNOTVARPSBC10["Annotated variant sites\n(SBC10-private)"]
        HIGHGENES["Genes with\nHIGH-impact variants"]
        OTHERGENES["Genes with\nnon-HIGH variants"]
  end
 subgraph s2["Epigenomics approach"]
        SBC10DMR["Differential methylation\nSBC10 vs others"]
        DMRGENES["Genes with\nDMR-ed promoter"]
  end
    VARPSBC10 -- SnpEFF annotation --> ANNOTVARPSBC10
    ANNOTVARPSBC10 -- "HIGH-impact filter" --> HIGHGENES
    ANNOTVARPSBC10 -- Other variants --> OTHERGENES
    SBC10DMR -- "Promoter annotation\n(2 kbp upstream, strand-aware)" --> DMRGENES
    HIGHGENES --> TIER1["Tier 1\nHIGH-impact variant + DMR"] & TIER2A["Tier 2\nHIGH-impact variant only"]
    DMRGENES --> TIER1 & TIER2B["Tier 2\nDMR only"] & TIER3["Tier 3\nnon-HIGH variant + DMR"]
    OTHERGENES --> TIER3 & TIER4["Tier 4\nnon-HIGH variant only"]

    style TIER1  fill:#5594dc
    style TIER2A fill:#BBDEFB
    style TIER2B fill:#BBDEFB
    style TIER3  fill:#dff0ff
    style TIER4  fill:#f5f5f5
```

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

3.  `annotate_vcf.sh`: Annotate private variants using SnpEFF. <br>
    ```shell
    ./docker/run.sh analysis/scripts/annotate_vcf.sh SBC10.private analysis/data/vcf/private_variants/SBC10.private.vcf.gz analysis/data/vcf/private_variants
    ```

4.  `annot_single_vcf_to_tsv.py`: Parse single-sampled VCF into TSV to explore in a notebook environment. The other file is designed to parse multi-sample (4 samples merged) VCF. <br>
    ```shell
    ./docker/run.sh python3 analysis/scripts/annot_single_vcf_to_tsv.py -v analysis/data/vcf/private_variants/SBC10.private.annotated.vcf.gz -o analysis/data/tsv
    ```

5.  `03_TAA/TAA.ipynb`: Analyze in notebook.

## Comparative epigenomics analysis
Technical steps
1.  `prepare_taa_dss.py`: Convert whole-genome filtered bedMethyl files to DSS input format. Collapses Watson/Crick CpG strand pairs for 5mC. <br>
    ```shell
    ./docker/run.sh python3 analysis/scripts/prepare_taa_dss.py
    ```
    Output: `analysis/data/taa_DSS/{sample}.5mC.dss.txt`

2.  `run_dss_dmr_taa.R`: Run DSS DML test and DMR calling for SBC10 vs each other accession. SBC10 is fixed as group1 so `diff.Methy > 0` means hypermethylated in SBC10. <br>
    ```shell
    ./docker/run.sh Rscript analysis/scripts/run_dss_dmr_taa.R
    ```
    Output: `analysis/data/taa_DMR/{pair}.5mC.DML.tsv`, `{pair}.5mC.DMR.tsv`, `DMR_summary.tsv`

3.  `annotate_DMR.py`: Annotate DMRs with overlapping genes, using strand-aware 2 kbp upstream promoter regions. <br>
    ```shell
    ./docker/run.sh python3 analysis/scripts/annotate_DMR.py
    ```
    Output: `analysis/data/taa_DMR/{pair}.5mC.DMR.annotated.tsv`

4.  `03_TAA/TAA.ipynb`: Integrate DMR gene lists with private variant gene lists to classify candidates by tier.