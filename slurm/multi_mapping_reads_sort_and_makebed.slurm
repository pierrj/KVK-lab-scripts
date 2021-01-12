#!/bin/bash
#SBATCH --job-name=multi_mapping_reads_G3_1A
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_with_mapq0/G3_1A_multimapping
SAMPLE=G3_1A

samtools sort -n ${SAMPLE}.sorted.mergedandpe.bwamem.bam -o name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam
bedtools bamtobed -cigar -i name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam > name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bed