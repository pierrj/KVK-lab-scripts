#!/bin/bash
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
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

## USAGE ##
# this script aligns two genomes using nucmer, then uses a python script to parse alignments that look like ecc mediated svs
# -r reference genome
# -q query genome
# -t threads
# -n prefix for output files
# -p filter for mummer plots based off alignment coverage percentage
# -o directory of output files
# -m temporary directory

ref=$(basename ${REFERENCE})
quer=$(basename ${QUERY})

if [ -d "${TMP_DIR}" ]; then
    rm -r ${TMP_DIR}
fi
mkdir ${TMP_DIR}

# get genome size files for query and reference
if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx ${REFERENCE}
fi
cut -f1,2 ${REFERENCE}.fai > ${TMP_DIR}/${ref}.genomesize
if [ ! -f "${QUERY}.fai" ]; then
    samtools faidx ${QUERY}
fi
cut -f1,2 ${QUERY}.fai > ${TMP_DIR}/${quer}.genomesize
# nucmer align genomes, output coords file for parsing
nucmer -t ${THREADS} --maxmatch -p ${TMP_DIR}/${OUTPUT_NAME} ${REFERENCE} ${QUERY}
show-coords ${TMP_DIR}/${OUTPUT_NAME}.delta > ${TMP_DIR}/${OUTPUT_NAME}.coords
tail -n+6 ${TMP_DIR}/${OUTPUT_NAME}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${TMP_DIR}/${OUTPUT_NAME}.processed_output.coords

# parse coords, output svs that look like they may have been caused by eccdnas
python /global/home/users/pierrj/git/python/ecc_sv_finder_AOC.py ${TMP_DIR}/${OUTPUT_NAME}.processed_output.coords ${TMP_DIR}/${ref}.genomesize ${TMP_DIR}/${quer}.genomesize ${TMP_DIR} ${ref} ${quer}

# generate jpgs of alignments and organize output files
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