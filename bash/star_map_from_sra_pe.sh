#!/bin/bash
while getopts g:1:2:s:a:b:t: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
2) READTWO=${OPTARG};;
s) SAMPLE=${OPTARG};;
a) ADAPTER_ONE=${OPTARG};;
b) ADAPTER_TWO=${OPTARG};;
t) THREADS=${OPTARG};;
esac
done

/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SAMPLE}.sra
cutadapt -a ${ADAPTER_ONE} -A ${ADAPTER_TWO} -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}
STAR --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DB} \
    --readFilesIn tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq