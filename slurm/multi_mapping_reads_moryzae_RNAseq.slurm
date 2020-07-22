#!/bin/bash
#SBATCH --job-name=multi_mapping_reads_moryzae_RNAseq
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
cd /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes

READONE='SRR298849.sra.fastq'
THREADS='24'
SAMPLE='SRR298849'
GENOME_DB='/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa'
bwa mem -a -t ${THREADS} ${GENOME_DB} ${READONE} -o tmp.merged.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.merged.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam