#!/bin/bash
while getopts g:1:s:t:m: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
m) MAPFILE=${OPTARG};;
esac
done

cutadapt -j ${THREADS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o tmp.trimmed.seqprep.${SAMPLE}.fastq ${READONE}
bwa mem -q -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam
samtools view -b -F 256 filtered.sorted.${SAMPLE}.bam > no_secondary.filtered.sorted.${SAMPLE}.bam
samtools view -b -q 1 no_secondary.filtered.sorted.${SAMPLE}.bam > uniq.filtered.sorted.${SAMPLE}.bam

# COMMENTS MISSING
samtools sort -n -@ ${THREADS} filtered.sorted.${SAMPLE}.bam > multimapped.filtered.name_sorted.${SAMPLE}.bam

rm tmp*
rm ${SAMPLE}.mergedandpe.bwamem.bam
rm filtered.sorted.${SAMPLE}.bam
rm ${SAMPLE}.sorted.mergedandpe.bwamem.bam