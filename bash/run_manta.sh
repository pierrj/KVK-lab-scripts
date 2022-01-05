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

/global/scratch/users/pierrj/manta-1.6.0.centos6_x86_64/bin/configManta.py \
    --runDir ${OUTDIR}/manta_out \
    --reference ${REFERENCE} \
    --bam ${FILE}

${OUTDIR}/manta_out/runWorkflow.py \
    --quiet \
    -m local \
    -j ${THREADS}

zcat ${OUTDIR}/manta_out/results/variants/candidateSV.vcf.gz > ${OUTPUT}