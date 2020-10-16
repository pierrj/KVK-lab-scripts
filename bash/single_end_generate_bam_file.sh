#!/bin/bash
while getopts g:1:s:t:m:q: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
m) MAPFILE=${OPTARG};;
q) READONE_LINK=${OPTARG};;
esac
done

curl -O ${READONE_LINK}
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}.fastq ${READONE}
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

rm tmp*
rm ${SAMPLE}.mergedandpe.bwamem.bam