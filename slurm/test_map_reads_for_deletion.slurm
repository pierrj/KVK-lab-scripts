#!/bin/bash
#SBATCH --job-name=test_map_reads_for_deletion
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/deletion_troubleshooting_test_map_reads

GENOME=13FM-16-1
FASTA=/global/scratch/users/pierrj/fungap_runs/gladieux_all/genomes_to_annotate/${GENOME}_genomic.fna
ACCESSION=SRR6669189
READONE=/global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${ACCESSION}/data/${ACCESSION}.sra_1.fastq
READTWO=/global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${ACCESSION}/data/${ACCESSION}.sra_2.fastq

# bwa index $FASTA

bwa mem -t ${SLURM_NTASKS} ${FASTA} ${READONE} ${READTWO} -o ${ACCESSION}_mapped_${GENOME}.bam

samtools sort ${ACCESSION}_mapped_${GENOME}.bam -o ${ACCESSION}_mapped_${GENOME}.sorted.bam

samtools index ${ACCESSION}_mapped_${GENOME}.sorted.bam