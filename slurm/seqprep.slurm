#!/bin/bash
#SBATCH --job-name=ecc_caller_G3_3A_seqprep
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/pureculture_samples/G3_3A/G3_3A
sample=G3_3A
/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f trimmed.${sample}_R1.fastq -r trimmed.${sample}_R2.fastq -1 seqprep.trimmed.${sample}_R1.fastq -2 seqprep.trimmed.${sample}_R2.fastq -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s merged.trimmed.${sample}