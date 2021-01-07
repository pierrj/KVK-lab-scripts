#!/bin/bash
#SBATCH --job-name=multi_mapping_reads_G3_1A
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_with_mapq0/G3_1A_multimapping

READONE='/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R1.fastq'
READTWO='/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R2.fastq'
THREADS='24'
SAMPLE='G3_1A'
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}
bwa mem -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam

# sort and index for filtering
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# filter using names of contigs in mapfile
samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

# remove tmp files and unfiltered files
# tmp command should probably be more specific, or use a tmp directory
rm tmp*
rm ${SAMPLE}.mergedandpe.bwamem.bam