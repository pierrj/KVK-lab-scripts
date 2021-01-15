#!/bin/bash
#SBATCH --job-name=sort_tests
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/samtools_sort_speed_tests

SAMPLE=G3_1A
READONE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R1.fastq
READTWO=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R2.fastq
THREADS=24
GENOME_DB=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}

# map reads
bwa mem -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam

# generate bam file from sam, this should probably just be passed as an option to bwa mem
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam

echo 'regular sort by coordinate'
date +"%T"

samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"

rm ${SAMPLE}.sorted.mergedandpe.bwamem.bam


echo 'multi threaded sort by coordinate'
date +"%T"

samtools sort -@ ${THREADS} ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"


echo 'multi threaded sort by qname'
date +"%T"

samtools sort -@ ${THREADS} ${SAMPLE}.mergedandpe.bwamem.bam > qname.${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"

rm qname.${SAMPLE}.sorted.mergedandpe.bwamem.bam

echo 'multi threaded sort by qname from coordinate sorted'
date +"%T"

samtools sort -@ ${THREADS} ${SAMPLE}.sorted.mergedandpe.bwamem.bam > qname.${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"