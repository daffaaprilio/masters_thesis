# Guide on running Clair3 in GPU mode

## Preparing Clair3 in `okadama`
> References <br> 
https://github.com/HKU-BAL/Clair3#installation, Option 1. Docker GPU (NVIDIA CUDA on Linux) is preferred <br>
https://docs.docker.com/get-started/docker_cheatsheet.pdf for list of useful Docker/Podman commands.

### Install libraries needed for accessing GPU
```shell
sudo dnf install nvidia-container-toolkit
sudo nvidia-ctk cdi generate --output=/etc/cdi/nvidia.yaml

# pull Clair3 image
podman pull docker.io/hkubal/clair3:v2.0.0_gpu
# check if image is pulled
podman images
```
### Trial run of the pulled image
```shell
# run image. Trial with just running bash
podman run -it \
  --device /dev/nvidia0 \
  --device /dev/nvidiactl \
  --device /dev/nvidia-uvm \
  --device /dev/nvidia-uvm-tools \
  -e CUDA_VISIBLE_DEVICES=0,1,2,3 \
  -v /home/daffa/Work/2026/thesis/resources/:/home/daffa/Work/2026/thesis/resources/ \
  -v /home/daffa/Work/2026/thesis/results/variant_calling/SBC10/:/home/daffa/Work/2026/thesis/results/variant_calling/SBC10/ \
  hkubal/clair3:v2.0.0_gpu \
  bash
```
### Actual running
```shell
# Running inside a screen session is recommended
podman run -it \
  --device /dev/nvidia0 \
  --device /dev/nvidiactl \
  --device /dev/nvidia-uvm \
  --device /dev/nvidia-uvm-tools \
  -e CUDA_VISIBLE_DEVICES=0,1,2,3 \
  -v /home/daffa/Work/2026/thesis/resources/:/home/daffa/Work/2026/thesis/resources/ \
  -v /home/daffa/Work/2026/thesis/results/variant_calling/SBC10/:/home/daffa/Work/2026/thesis/results/variant_calling/SBC10/ \
  hkubal/clair3:v2.0.0_gpu \
  /opt/bin/run_clair3.sh \
    --bam_fn=/home/daffa/Work/2026/thesis/resources/align_bam_sample/SBC10.bam \
    --ref_fn=/home/daffa/Work/2026/thesis/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
    --threads=8 \
    --platform=ont \
    --model_path=/opt/models/r1041_e82_400bps_sup_v520_with_mv \
    --output=/home/daffa/Work/2026/thesis/results/variant_calling/SBC10/ \
    --use_gpu \
    --include_all_ctgs
```
### Alternative approach
Building dockerfile consists of the Clair3 image and includes necessary library, `nvidia-container-toolkit`
```shell
podman build -f Dockerfile.clair3 -t clair3-thesis:latest .

podman run -it --rm \
  --device nvidia.com/gpu=all \
  -e CUDA_VISIBLE_DEVICES=0,1,2,3 \
  -v /home/daffa/Work/2026/thesis/resources/:/home/daffa/Work/2026/thesis/resources/ \
  -v /home/daffa/Work/2026/thesis/results/:/home/daffa/Work/2026/thesis/results/ \
  clair3-thesis:latest \
  /opt/bin/run_clair3.sh \
    --bam_fn=/home/daffa/Work/2026/thesis/resources/align_bam_sample/SBC10.bam \
    --ref_fn=/home/daffa/Work/2026/thesis/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
    --threads=8 --platform=ont \
    --model_path=/opt/models/r1041_e82_400bps_sup_v520_with_mv \
    --output=/home/daffa/Work/2026/thesis/results/variant_calling/SBC10/ \
    --use_gpu --include_all_ctgs
```