# Multi-omics Integration for Sorghum Metabolic Variation Analysis
Consists of three omics analyses: gene co-expression analysis (transcriptomics), variant analysis (genomics) and methylation analysis (epigenomics).

## Preparation
Install the following tools:
- Snakemake
- Samtools
- BCFtools
- HTSlib

## Genomics
### Read data preparation
Jobs included:
- Read alignment
- Aligned reads indexing
- Reads coverage/depth calculation <br>

All these jobs are handled/included in `workflow/rules/reads_preprocessing.smk`. <br>

`reads_preprocessing.smk` processes each read as a library. There are 6 libraries:
| Library | Sample | TAA Production |
|---------|--------|-----------|
| r0074 | SBC4 | ++  |
| r0066 | SBC10 | +++ |
| r0075 <br> r0078 <br> r0078-2 | SBC11 | - |
| r0076 | SBC23 | ++ |

Running `reads_preprocessing.smk` Snakefile
```shell
# Running from the beginning
snakemake -s workflow/rules/reads_preprocessing.smk -c 24 -j 6 -pn
# Picking up from the pickle file, 
snakemake -s workflow/rules/reads_preprocessing.smk just_plot -c 24 -j 6 -pn
```

#### Converting library to sample
For SBC4, SBC10, and SBC23, use symbolic link to save storage, preventing file duplication.
```shell
declare -A samples=(
    [SBC10]="r0066"
    [SBC4]="r0074"
    [SBC23]="r0076"
)
for sample in ${(k)samples}; do
    ln -s "resources/align_bam/${samples[$sample]}.bam" "resources/align_bam_sample/${sample}.bam"
done
```
For SBC11 libraries:
```shell
# merge multiple libraries
samtools merge -r -@ 6 resources/align_bam_sample/SBC11.bam resources/align_bam/r0075.bam resources/align_bam/r0078.bam resources/align_bam/r0078-2.bam
# sort and index merged bam
samtools sort -@ 6 -o resources/align_bam_sample/SBC11.bam resources/align_bam_sample/SBC11.bam
samtools index resources/align_bam_sample/SBC11.bam
```
Continue making depth file and its statistics for SBC11
```shell
# obtain samtools .depth file
samtools depth -a resources/align_bam_sample/SBC11.bam -o resources/depth/SBC11.depth
# plot the depth file — saves SBC11_depth.png, .svg, and .pdf
python workflow/scripts/visualize_depth.py \
    --input resources/depth/SBC11.depth \
    --output resources/depth/SBC11_depth \
    --library "r0078, r0075, r0078-2" \
    --save-pickle resources/depth/SBC11_depth.pkl
# re-plot from existing pickle — saves SBC11_depth.png, .svg, and .pdf
python workflow/scripts/visualize_depth.py \
    --load-pickle resources/depth/SBC11_depth.pkl \
    --output resources/depth/SBC11_depth \
    --library "r0078, r0075, r0078-2"
```

### Variant calling with naive model
Clair3 is used for variant calling. For the first stage, variant calling is done without pre-training the model on sorghum data. Pre-training the model with sorghum data is also considered.

Variant calling using Snakefile
```shell
# designed to run clair3_cpu rule
snakemake -s workflow/rules/variant_analysis.smk -c 24 -j 4 -pn
# same, but only on 2 top samples
snakemake --snakefile workflow/rules/variant_analysis.smk \
  /home/daffa/Work/2026/thesis/results/variant_calling/SBC10/merge_output.vcf.gz \
  /home/daffa/Work/2026/thesis/results/variant_calling/SBC11/merge_output.vcf.gz \
  -c 8 -j 2 -pn
```

Tidy up `Clair3`'s output directory
```shell
# symbolic link to results/vcf/SBC10.vcf.gz
cd /home/daffa/Work/2026/thesis/results/vcf
ln -s ../variant_calling/SBC10/merge_output.vcf.gz SBC10.vcf.gz
# index VCF
bcftools index SBC10.vcf.gz 
```

### VCF postprocessing
Download SnpEFF database for Sorghum
```shell
snpEff download Sorghum_bicolor
```
Downloaded database is saved in `/home/daffa/local/bin/snpEff/data/Sorghum_bicolor/` <br>
Prepare the chromosome conversion in the database. The trick: create a conversion key by reading the following files:
- any VCF file, we need the header
```shell
bcftools view -h results/vcf/SBC4.vcf.gz | less -S
```
This has been done and saved in `workflow/scripts/vcf_chr_list.txt`
```
##contig=<ID=NC_012870.2,length=80884392>
##contig=<ID=NC_012871.2,length=77742459>
##contig=<ID=NC_012872.2,length=74386277>
##contig=<ID=NC_012873.2,length=68658214>
```
- the SnpEff database
```shell
snpEff dump Sorghum_bicolor | less -S
```
This, also has been done and saved in `workflow/scripts/snpeff_db_chr_list.txt`, but first, don't forget to trim the unnecessary part (from the last contigs to the end of the file) to prevent the file from bloating.
```
#-----------------------------------------------
# Number of chromosomes      : 867
# Chromosomes                : Format 'chromo_name size codon_table'
#		'1'	80884392	Standard
#		'2'	77742459	Standard
#		'3'	74386277	Standard
The script below will automatically create the conversion key
```shell
python3 workflow/scripts/creating_synonyms_chr.py
```
Run Snakefile for VCF file processing
```shell
snakemake --snakefile workflow/rules/vcf_processing.smk -c 24 -j 4 -pn
```
Overview of rules in this file:
- Adding SAMPLE information to each VCF　(`bcftools reheader`)
- Filter VCF, only pass variants from 10 chromosomes with high quality
- Phase (add haplotype information) to filtered VCF
- Performs consequences calling (`snpEff ann`)
