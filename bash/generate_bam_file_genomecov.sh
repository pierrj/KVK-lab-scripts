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

/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${READONE} -r ${READTWO} -1 tmp.seqprep.${READONE} -2 tmp.seqprep.${READTWO} -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${READONE} -p tmp.trimmed.seqprep.${READTWO} tmp.seqprep.${READONE} tmp.seqprep.${READTWO}
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --nextseq-trim=20 -o tmp.trimmed.merged tmp.merged
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${READONE} tmp.trimmed.seqprep.${READTWO} -o tmp.seqprep.trimmed.bwamem.sam
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.merged -o tmp.merged.trimmed.bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.bwamem.sam > tmp.seqprep.trimmed.bwamem.bam
samtools view -S -b tmp.merged.trimmed.bwamem.sam > tmp.merged.trimmed.bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.bwamem.bam -in tmp.merged.trimmed.bwamem.bam -out ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

bedtools genomecov -d -ibam ${SAMPLE}.sorted.mergedandpe.bwamem.bam > ${SAMPLE}.genomecoverage.bed

rm tmp*