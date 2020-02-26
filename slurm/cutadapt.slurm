#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ../tests/trimmed.IF_1A_S76_L002_R1_001.fastq.gz -p ../tests/IF_1A_S76_L002_R2_001.fastq.gz IF_1A_S76_L002_R1_001.fastq.gz IF_1A_S76_L002_R2_001.fastq.gz