#!/bin/bash
#SBATCH --job-name=filteringbymapping
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
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_3B
samtools view -b -h -f 4 pBR322_CT_3B.sam > pBR322_filtered_CT_3B.bam
samtools sort -n pBR322_filtered_CT_3B.bam > sorted.pBR322_filtered_CT_3B.bam
bedtools bamtofastq -i sorted.pBR322_filtered_CT_3B.bam -fq sorted.pBR322_filtered_CT_3B_R1.fastq -fq2 sorted.pBR322_filtered_CT_3B_R2.fastq
fastqc sorted.pBR322_filtered_CT_3B_R1.fastq
fastqc sorted.pBR322_filtered_CT_3B_R2.fastq