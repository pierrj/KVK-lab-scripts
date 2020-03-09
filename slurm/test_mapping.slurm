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
cd //global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/pureculture_samples/G3_1A
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/pBR322.fasta G3_1A_R1.fastq G3_1A_R2.fastq -o pBR322_G3_1A.sam
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/stress_experiments/references/B71_mitochondria.fasta G3_1A_R1.fastq G3_1A_R2.fastq -o B71_G3_1A.sam