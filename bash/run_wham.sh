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


/global/scratch/users/pierrj/wham/bin/whamg -f ${FILE} -a ${REFERENCE} -x ${THREADS} > ${OUTPUT}