# Multi-omics Integration for Sorghum Metabolic Variation Analysis
Consists of three omics analyses: gene co-expression analysis (transcriptomics), variant analysis (genomics) and methylation analysis (epigenomics).

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

### Running `reads_preprocessing.smk` Snakefile
```shell
# Running from the beginning
snakemake -s workflow/rules/reads_preprocessing.smk -c 24 -j 6 -pn
# Picking up from the pickle file, 
snakemake -s workflow/rules/reads_preprocessing.smk just_plot -c 24 -j 6 -pn
```

### Converting library to sample
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
# plot the depth file
python workflow/scripts/visualize_depth.py \
    --input resources/depth/SBC11.depth \
    --output resources/depth/SBC11_depth.svg \
    --library "r0078, r0075, r0078-2" \
    --save-pickle resources/depth/SBC11_depth.pkl
# plot the depth file (by reading existing SBC11_depth.pkl file)
python workflow/scripts/visualize_depth.py \
    --load-pickle resources/depth/SBC11_depth.pkl \
    --output resources/depth/SBC11_depth.svg \
    --library "r0078, r0075, r0078-2" 
```

### Variant calling with naive model
Clair3 is used for variant calling. For the first stage, variant calling is done without pre-training the model on sorghum data. Pre-training the model with sorghum data is also considered.

#### Transferring files from `matsu` to `okadama`
```shell
# reference genome
rsync --dry-run -av -s daffa@matsu:/home/daffa/Work/2026/thesis/resources/ref/* resources/ref/
# aligned bam (both sample level and library)
rsync --dry-run -av -s daffa@matsu:/home/daffa/Work/2026/thesis/resources/align_bam/* resources/align_bam/
rsync --dry-run -av -s daffa@matsu:/home/daffa/Work/2026/thesis/resources/align_bam_sample/* resources/align_bam_sample/
```

#### Preparing Clair3 in `okadama`
> References <br> 
https://github.com/HKU-BAL/Clair3#installation, Option 1. Docker GPU (NVIDIA CUDA on Linux) is preferred <br>
https://docs.docker.com/get-started/docker_cheatsheet.pdf for list of useful Docker/Podman commands.
```shell
# pull Clair3 image
podman pull docker.io/hkubal/clair3:v2.0.0_gpu
# check if image is pulled
podman images
# run image. Running inside a screen session is recommended
podman run -it \
  --device /dev/nvidia0 \
  --device /dev/nvidiactl \
  --device /dev/nvidia-uvm \
  --device /dev/nvidia-uvm-tools \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:v2.0.0_gpu \
  /opt/bin/run_clair3.sh \
    --bam_fn=${INPUT_DIR}/input.bam \
    --ref_fn=${INPUT_DIR}/ref.fa \
    --threads=${THREADS} \
    --platform=ont \
    --model_path=/opt/models/${MODEL_NAME} \
    --output=${OUTPUT_DIR} \
    --use_gpu
```

### Draft genome assembly
