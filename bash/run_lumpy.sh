#!/bin/bash
while getopts f:r:t:o:d: option
do
case "${option}"
in
f) FILE=${OPTARG};;
r) REFERENCE=${OPTARG};;
t) THREADS=${OPTARG};;
o) OUTPUT=${OPTARG};;
d) OUTDIR=${OPTARG};;
esac
done

basename_file=$(basename ${FILE} )

smoove call --name ${basename_file} \
    --fasta ${REFERENCE} \
    --processes ${THREADS} \
    --outdir ${FILE}_smoove_out \
    ${FILE}

zcat ${FILE}_smoove_out/${basename_file}-smoove.vcf.gz > ${OUTPUT}