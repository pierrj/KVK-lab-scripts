#!/bin/bash
#SBATCH --job-name=test_sv_calling_workshop
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv_calling_workshop

module load sra-tools

prefetch ERR3624573 -O .
fasterq-dump -e 20 -O . -t tmp ERR3624573/ERR3624573.sra

ACCESSION=ERR3624573
REFERENCE_GENOME=GCA_000001735.2_TAIR10.1_genomic.fna

module load bwa

bwa index ${REFERENCE_GENOME}

bwa mem -R "@RG\tID:${ACCESSION}\tSM:${ACCESSION}" \
    -t 24 \
    ${REFERENCE_GENOME} \
    ${ACCESSION}_1.fastq \
    ${ACCESSION}_2.fastq \
    -o ${ACCESSION}.bam

source activate /global/scratch/users/pierrj/conda_envs/lumpy

smoove call --name ERR3624573 \
	--fasta GCA_000001735.2_TAIR10.1_genomic.fna \
	--processes 24 \
	--outdir ERR3624573_smoove_out \
	ERR3624573.bam

zcat ERR3624573_smoove_out/ERR3624573-smoove.vcf.gz > ERR3624573.lumpy.vcf