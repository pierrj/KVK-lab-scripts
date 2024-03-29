#!/bin/bash
while getopts g:1:s:t: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
esac
done

/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SAMPLE}.sra
STAR --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DB} \
    --readFilesIn ${READONE} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${SAMPLE}. \
    --outSAMstrandField intronMotif
cufflinks -p ${THREADS} --library-type ff-unstranded -o cufflinks.${sample}.dir ${SAMPLE}.Aligned.sortedByCoord.out.bam