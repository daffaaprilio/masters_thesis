# Multi-omics integration for sorghum metabolic variation
Consists of 3 omics analysis: gene co-expression analysis (transcriptomics), variant analysis (genomics) and methylation analysis (epigenomics).
## Genomics
Preparing read files for variant calling.
```shell
# run snakemake (current version only designed for variant calling)
snakemake -s workflow/rules/reads_preprocessing.smk -c 24 -j 6 -np
```
