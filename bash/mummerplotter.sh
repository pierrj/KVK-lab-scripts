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
while getopts r:q:e:u:o:p: option
do
case "${option}"
in
r) REFERENCE=${OPTARG};;
q) QUERY=${OPTARG};;
e) ref=${OPTARG};;
u) quer=${OPTARG};;
o) SUBDIR_OUTPUT=${OPTARG};;
p) PERCENT_ZEROES_FILTER=${OPTARG};;
esac
done

## USAGE ##
# this script generates mummer plots between two bed files and genome files
# -r reference fasta file
# -q query fasta file
# -e ref prefix
# -u quer prefix
# -o output subdirectory
# -p percent of alignments not covered by sequences

output_realpath=$(realpath ${SUBDIR_OUTPUT})
SUBDIR_TMP=tmp
if [ -d "${SUBDIR_TMP}" ]; then
    rm -r ${SUBDIR_TMP}
fi
mkdir ${SUBDIR_TMP}
# number of alignments to process
alignment_number=$(ls -1 | wc -l | awk '{print ($1-1)/2-1}')
seq 0 ${alignment_number} > ${SUBDIR_TMP}/mapfile
while read plot; do
    # grab bed files and get fasta sequences from them
    bedtools getfasta -fi ${REFERENCE} -bed ${ref}_${plot}.bed -fo ${SUBDIR_TMP}/${ref}_${plot}.fasta
    samtools faidx ${SUBDIR_TMP}/${ref}_${plot}.fasta
    # get genome size files
    cut -f1,2 ${SUBDIR_TMP}/${ref}_${plot}.fasta.fai > ${SUBDIR_TMP}/${ref}_${plot}.genomesize
    bedtools getfasta -fi ${QUERY} -bed ${quer}_${plot}.bed -fo ${SUBDIR_TMP}/${quer}_${plot}.fasta
    # align fastas
    nucmer -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${ref}_${plot}.fasta ${SUBDIR_TMP}/${quer}_${plot}.fasta
    show-coords ${SUBDIR_TMP}/${plot}.delta | tail -n +6 | awk -v OFS='\t' '{print $12, $1, $2}' | sort -k1,1 -k2,2n >  ${SUBDIR_TMP}/${plot}.bed
    # count percentage of alignment not covered
    bedtools genomecov -d -i ${SUBDIR_TMP}/${plot}.bed -g ${SUBDIR_TMP}/${ref}_${plot}.genomesize > ${SUBDIR_TMP}/${plot}.genomecov
    total_size=$(wc -l ${SUBDIR_TMP}/${plot}.genomecov | awk '{print $1}')
    size_zeroes=$(awk '$3==0' ${SUBDIR_TMP}/${plot}.genomecov | wc -l | awk '{print $1}')
    percent_zeroes=$(awk -v var1=$size_zeroes -v var2=$total_size 'BEGIN { OFMT="%f";print  ( var1 / var2 ) }')
    if (( $(echo "$percent_zeroes < ${PERCENT_ZEROES_FILTER}" |bc -l) )) # filter based off percentage
    then
        # draw mummer plot, make pdf and convert to jpg
        mummerplot --color -postscript -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${plot}.delta
        gnuplot ${SUBDIR_TMP}/${plot}.gp
        ps2pdf ${SUBDIR_TMP}/${plot}.ps ${SUBDIR_TMP}/${plot}.pdf
        convert -density 150 ${SUBDIR_TMP}/${plot}.pdf -quality 90 ${output_realpath}/${ref}_v_${quer}_${plot}.jpg
        cp ${SUBDIR_TMP}/${plot}.pdf ${output_realpath}/${ref}_v_${quer}_${plot}.pdf
    fi
done < ${SUBDIR_TMP}/mapfile