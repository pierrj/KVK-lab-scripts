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

delly call \
    -g ${REFERENCE} \
    -o ${OUTDIR}/delly.out \
    ${FILE}

bcftools filter \
    -O v \
    -o ${OUTPUT} \
    -i "FILTER == 'PASS'" \
    ${OUTDIR}/delly.out