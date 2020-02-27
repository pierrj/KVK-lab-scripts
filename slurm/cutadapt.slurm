#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed.CT_3C_S57_L002_R1_001.fastq.gz -p trimmed.CT_3C_S57_L002_R2_001.fastq.gz CT_3C_S57_L002_R1_001.fastq.gz CT_3C_S57_L002_R2_001.fastq.gz
fastqc trimmed.CT_3C_S57_L002_R1_001.fastq.gz --outdir=fastqc_output