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
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/B71_mitochondria.fasta CT_3B_R1.fastq CT_3B_R2.fastq -o B71_CT_3B.sam
samtools view -b B71_CT_3B.sam > B71_CT_3B.bam
samtools sort B71_CT_3B.bam -o sorted.B71_CT_3B.bam