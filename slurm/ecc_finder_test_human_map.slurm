#!/bin/bash
#SBATCH --job-name=ecc_finder_test_human_map
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


SAMPLE=SRR6315399
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"

cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${SAMPLE}/ecc_finder


source activate /global/scratch/users/pierrj/conda_envs/ecc_finder_pre_steps

module purge


fastp --thread ${SLURM_NTASKS} -i /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${SAMPLE}/${SAMPLE}.sra_1.fastq \
            -I /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${SAMPLE}/${SAMPLE}.sra_2.fastq \
            -o ecc.sr.R1.fastq.gz \
            -O ecc.sr.R2.fastq.gz \
            -h ecc.sr.html

bwa mem -t ${SLURM_NTASKS} ${genome_bwa} \
            ecc.sr.R1.fastq.gz \
            ecc.sr.R2.fastq.gz > ecc.sr.sam 2> ecc.sr.sam.log