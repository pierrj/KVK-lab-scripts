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

prefetch ERR3415817 -O .
fasterq-dump -e 20 -O . -t tmp ERR3415817.sra

ACCESSION=ERR3415817
REFERENCE_GENOME=GCA_000001735.2_TAIR10.1_genomic.fna

module load bwa

bwa index ${REFERNCE_GENOME}

bwa mem -R "@RG\tID:${ACCESSION}\tSM:${ACCESSION}" \ # read group labels
	-t 24 \ # number of threads, should be same as ntasks in header
	${REFERENCE_GENOME} \ # reference
	${ACCESSION}.sra_1.fastq \ # read file one
	${ACCESSION}.sra_2.fastq \ # read file two
	-o ${ACCESSION}.bam # output mapped reads file

source activate /global/scratch/users/pierrj/conda_envs/lumpy

smoove call --name ERR3415817 \
	--fasta GCA_000001735.2_TAIR10.1_genomic.fna \
	--processes 24 \
	--outdir ERR3415817_smoove_out \
	ERR3415817.bam

zcat ERR3415817_smoove_out/ERR3415817-smoove.vcf.gz > ERR3415817.lumpy.vcf