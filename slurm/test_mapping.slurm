#!/bin/bash
#SBATCH --job-name=test_mapping
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
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/pBR322.fasta CT_3B_R1.fastq CT_3B_R2.fastq -o pBR322_CT_3B.sam
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/B71_mitochondria.fasta CT_3B_R1.fastq CT_3B_R2.fastq -o B71_CT_3B.sam
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/HP_3A
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/pBR322.fasta HP_3A_R1.fastq HP_3A_R2.fastq -o pBR322_HP_3A.sam
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/B71_mitochondria.fasta HP_3A_R1.fastq HP_3A_R2.fastq -o B71_HP_3A.sam