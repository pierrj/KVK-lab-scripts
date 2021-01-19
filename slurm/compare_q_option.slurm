#!/bin/bash
#SBATCH --job-name=compare_q_option
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

READONE='/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R1.fastq'
READTWO='/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R2.fastq'
THREADS='24'
SAMPLE='G3_1A'
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/compare_q_option/noqoption

bwa mem -a -t ${THREADS} ${GENOME_DB} ${READONE} ${READTWO} -o test_mapping.sam

samtools view -c test_mapping.sam

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/compare_q_option/yesqoption

bwa mem -q -a -t ${THREADS} ${GENOME_DB} ${READONE} ${READTWO} -o test_mapping.sam

samtools view -c test_mapping.sam