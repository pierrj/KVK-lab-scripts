#!/bin/bash
while getopts g:1:2:s:t: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
2) READTWO=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
esac
done

/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${READONE} -r ${READTWO} -1 tmp.seqprep.${SAMPLE}_R1.fastq.gz -2 tmp.seqprep.${SAMPLE}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${SAMPLE}.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --nextseq-trim=20 -o tmp.trimmed.merged.${SAMPLE}.fastq tmp.merged.${SAMPLE}.fastq.gz
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.merged.${SAMPLE}.fastq -o tmp.merged.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${SAMPLE}_bwamem.sam > tmp.merged.trimmed.${SAMPLE}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam -in tmp.merged.trimmed.${SAMPLE}_bwamem.bam -out ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

bedtools genomecov -d -ibam ${SAMPLE}.sorted.mergedandpe.bwamem.bam > ${SAMPLE}.genomecoverage.bed

rm tmp*