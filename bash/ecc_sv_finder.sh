#!/bin/bash
while getopts r:q:t:o:p: option
do
case "${option}"
in
r) REFERENCE=${OPTARG};;
q) QUERY=${OPTARG};;
t) THREADS=${OPTARG};;
n) OUTPUT_NAME=${OPTARG};;
p) PERCENT_ZEROES_FILTER=${OPTARG};;
o) OUTPUT_DIR=${OPTARG};;
esac
done

ref=$(basename ${REFERENCE})
quer=$(basename ${QUERY})

## two different versions of mummer used which is awk

samtools faidx ${REFERENCE}
cut -f1,2 ${REFERENCE}.fai > ${ref}.genomesize
samtools faidx ${QUERY}
cut -f1,2 ${QUERY}.fai > ${quer}.genomesize
/global/scratch/pierrj/mummer_4/bin/nucmer -t ${THREADS} --maxmatch -p ${OUTPUT_NAME} ${REFERENCE} ${QUERY}
/global/scratch/pierrj/mummer_4/bin/show-coords ${OUTPUT_NAME}.delta > ${OUTPUT_NAME}.coords
tail -n+6 ${OUTPUT_NAME}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${OUTPUT_NAME}.processed_output.coords

python /global/home/users/pierrj/git/python/ecc_sv_finder_AOC.py ${OUTPUT_NAME}.processed_output.coords ${ref}.genomesize ${quer}.genomesize . ${ref} ${quer}

## mummer plot stuff here

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir ${OUTPUT_DIR}
fi
cd ${ref}_v_${quer}
/global/home/users/pierrj/git/bash/mummerplotter.sh -r ${REFERENCE} -q ${QUERY}  -e ${ref} -u ${quer} -o ${OUTPUT_DIR} -p ${PERCENT_ZEROES_FILTER}