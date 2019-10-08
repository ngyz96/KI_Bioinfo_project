#!/bin/bash

#SBATCH -A snic2019-8-250
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 4:00:00
#SBATCH -J STAR

module unload
module load bioinfo-tools
module load star/2.7.2b

GENOME_FASTA='/domus/h1/ngyz96/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa'
GENOME_DIR='/proj/snic2019-8-250/RNOR_STAR_GENOME'
GENOME_GTF='/domus/h1/ngyz96/Rattus_norvegicus.Rnor_6.0.98.gtf'

STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${GENOME_DIR} --genomeFastaFiles ${GENOME_FASTA} --sjdbGTFfile ${GENOME_GTF} 
