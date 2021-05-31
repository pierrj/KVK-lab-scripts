#!/bin/bash
while getopts r:q:t:n:p:o:m: option
do
case "${option}"
in
r) REFERENCE=${OPTARG};;
q) QUERY=${OPTARG};;
t) THREADS=${OPTARG};;
n) OUTPUT_NAME=${OPTARG};;
p) PERCENT_ZEROES_FILTER=${OPTARG};;
o) OUTPUT_DIR=${OPTARG};;
m) TMP_DIR=${OPTARG};;
esac
done

ref=$(basename ${REFERENCE})
quer=$(basename ${QUERY})

## two different versions of mummer used which is awk

if [ -d "${TMP_DIR}" ]; then
    rm -r ${TMP_DIR}
fi
mkdir ${TMP_DIR}


if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx ${REFERENCE}
fi
cut -f1,2 ${REFERENCE}.fai > ${TMP_DIR}/${ref}.genomesize
if [ ! -f "${QUERY}.fai" ]; then
    samtools faidx ${QUERY}
fi
cut -f1,2 ${QUERY}.fai > ${TMP_DIR}/${quer}.genomesize
/global/scratch/pierrj/mummer_4/bin/nucmer -t ${THREADS} --maxmatch -p ${TMP_DIR}/${OUTPUT_NAME} ${REFERENCE} ${QUERY}
/global/scratch/pierrj/mummer_4/bin/show-coords ${TMP_DIR}/${OUTPUT_NAME}.delta > ${TMP_DIR}/${OUTPUT_NAME}.coords
tail -n+6 ${TMP_DIR}/${OUTPUT_NAME}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${TMP_DIR}/${OUTPUT_NAME}.processed_output.coords

python /global/home/users/pierrj/git/python/ecc_sv_finder_AOC.py ${TMP_DIR}/${OUTPUT_NAME}.processed_output.coords ${TMP_DIR}/${ref}.genomesize ${TMP_DIR}/${quer}.genomesize ${TMP_DIR} ${ref} ${quer}

## mummer plot stuff here

if [ -d "${TMP_DIR}/${ref}_v_${quer}" ]; then
    output_realpath=$(realpath ${OUTPUT_DIR})
    reference_realpath=$(realpath ${REFERENCE})
    query_realpath=$(realpath ${QUERY})
    if [ -d "${output_realpath}" ]; then
        rm -r ${output_realpath}
    fi
    mkdir ${output_realpath}
    cd ${TMP_DIR}/${ref}_v_${quer}
    /global/home/users/pierrj/git/bash/mummerplotter.sh -r ${reference_realpath} -q ${query_realpath}  -e ${ref} -u ${quer} -o ${output_realpath} -p ${PERCENT_ZEROES_FILTER}
fi