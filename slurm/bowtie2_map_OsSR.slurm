#!/bin/bash
#SBATCH --job-name=rice_bowtie2_mapping
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/OsSR_bowtie
SAMPLE=OsSR
READONE="OsSR_R1.fastq"
READTWO="OsSR_R2.fastq"

/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${READONE} -r ${READTWO} -1 tmp.seqprep.${SAMPLE}_R1.fastq.gz -2 tmp.seqprep.${SAMPLE}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${SAMPLE}.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --nextseq-trim=20 -o tmp.trimmed.merged.${SAMPLE}.fastq tmp.merged.${SAMPLE}.fastq.gz
bowtie2 --sentitive-local -p ${SLURM_NTASKS} -k 1 -x ORSA_w_organelles_bwt2 -1 tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -2 tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -S tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
bowtie2 --sentitive-local -p ${SLURM_NTASKS} -k 1 -x ORSA_w_organelles_bwt2 -U tmp.trimmed.merged.${SAMPLE}.fastq -S tmp.merged.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${SAMPLE}_bwamem.sam > tmp.merged.trimmed.${SAMPLE}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam -in tmp.merged.trimmed.${SAMPLE}_bwamem.bam -out ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

bedtools genomecov -d -ibam ${SAMPLE}.sorted.mergedandpe.bwamem.bam > ${SAMPLE}.genomecoverage.bed