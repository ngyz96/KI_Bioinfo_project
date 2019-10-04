#!/bin/bash -l

#SBATCH -A snic2019-8-250
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 20:00:00
#SBATCH -J fastqc

module unload
module load bioinfo-tools
module load FastQC/0.11.8
find /proj/uppstore2019102/ESCG_data/ -type f -name "P13104_1000?_S?_L00?_R?_001.fastq.gz" -print0 | xargs -0  -n 1 fastqc --outdir=/home/ngyz96/fastqc_result
