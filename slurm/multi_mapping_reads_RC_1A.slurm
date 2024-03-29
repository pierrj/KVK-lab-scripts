#!/bin/bash
#SBATCH --job-name=multi_mapping_reads_RC_1A
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/RC_1A_multimapping

READONE='/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R1.fastq'
READTWO='/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R2.fastq'
THREADS='24'
SAMPLE='RC_1A'
GENOME_DB='/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa'
/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${READONE} -r ${READTWO} -1 tmp.seqprep.${SAMPLE}_R1.fastq.gz -2 tmp.seqprep.${SAMPLE}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${SAMPLE}.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --nextseq-trim=20 -o tmp.trimmed.merged.${SAMPLE}.fastq tmp.merged.${SAMPLE}.fastq.gz
bwa mem -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
bwa mem -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.merged.${SAMPLE}.fastq -o tmp.merged.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${SAMPLE}_bwamem.sam > tmp.merged.trimmed.${SAMPLE}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam -in tmp.merged.trimmed.${SAMPLE}_bwamem.bam -out ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# bedtools genomecov -d -ibam ${SAMPLE}.sorted.mergedandpe.bwamem.bam > ${SAMPLE}.genomecoverage.bed