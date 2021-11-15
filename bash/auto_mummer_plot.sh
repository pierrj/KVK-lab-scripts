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
while getopts r:g:s:p: option
do
case "${option}"
in
r) REFERENCE=${OPTARG};;
g) GENOME_DIR=${OPTARG};;
s) SUBDIR_OUTPUT=${OPTARG};;
p) PERCENT_ZEROES_FILTER=${OPTARG};;
esac
done

## USAGE ##
# this script aligns a reference genome to a directory containing bed files of locations of regions in the reference and query genomes and aligns them to each other using nucmer
# it outputs the location of alignments as well as a jpeg of the alignment
# -r reference genome fasta
# -g location of query genomes directory
# -s name of output jpgs
# -p percentage of gaps of alignments to use to filter out outputs

if [ -d "${SUBDIR_OUTPUT}" ]; then
    rm -r ${SUBDIR_OUTPUT}
fi

ls -1 | grep -v genome_mapfile | awk '{print substr($1,9)}'> genome_mapfile

mkdir ${SUBDIR_OUTPUT}

while read genome ; do
    cd guy11_v_${genome}
    SUBDIR_TMP=tmp
    if [ -d "${SUBDIR_TMP}" ]; then
        rm -r ${SUBDIR_TMP}
    fi
    mkdir ${SUBDIR_TMP}
    alignment_number=$(ls -1 | wc -l | awk '{print ($1-1)/2-1}')
    seq 0 ${alignment_number} > ${SUBDIR_TMP}/mapfile
    ISOLATE=${GENOME_DIR}/${genome}_genomic.fna
    while read plot; do
        bedtools getfasta -fi ${REFERENCE} -bed guy11_${plot}.bed -fo ${SUBDIR_TMP}/guy11_${plot}.fasta
        samtools faidx ${SUBDIR_TMP}/guy11_${plot}.fasta
        cut -f1,2 ${SUBDIR_TMP}/guy11_${plot}.fasta.fai > ${SUBDIR_TMP}/guy11_${plot}.genomesize # genome size files necessary for downstream analysis....
        bedtools getfasta -fi ${ISOLATE} -bed ${genome}_${plot}.bed -fo ${SUBDIR_TMP}/${genome}_${plot}.fasta # get sequences to align from bedfiles
        nucmer -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/guy11_${plot}.fasta ${SUBDIR_TMP}/${genome}_${plot}.fasta # maxmatch is not included here so this is mostly for aligning short stretches of genomes to each other
        show-coords ${SUBDIR_TMP}/${plot}.delta | tail -n +6 | awk -v OFS='\t' '{print $12, $1, $2}' | sort -k1,1 -k2,2n >  ${SUBDIR_TMP}/${plot}.bed
        bedtools genomecov -d -i ${SUBDIR_TMP}/${plot}.bed -g ${SUBDIR_TMP}/guy11_${plot}.genomesize > ${SUBDIR_TMP}/${plot}.genomecov
        total_size=$(wc -l ${SUBDIR_TMP}/${plot}.genomecov | awk '{print $1}') # total size of alignment
        size_zeroes=$(awk '$3==0' ${SUBDIR_TMP}/${plot}.genomecov | wc -l | awk '{print $1}') # calculate number of gaps in alignment
        percent_zeroes=$(awk -v var1=$size_zeroes -v var2=$total_size 'BEGIN { OFMT="%f";print  ( var1 / var2 ) }') # percentage
        if (( $(echo "$percent_zeroes < ${PERCENT_ZEROES_FILTER}" |bc -l) )) # make sure percentage of gaps are lower than cutoff
        then
            mummerplot --color -postscript -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${plot}.delta # draw mummer plot and convert to jpg
            ps2pdf ${SUBDIR_TMP}/${plot}.ps ${SUBDIR_TMP}/${plot}.pdf
            convert -density 150 ${SUBDIR_TMP}/${plot}.pdf -quality 90 ../${SUBDIR_OUTPUT}/${genome}_${plot}.jpg
        fi
    done < ${SUBDIR_TMP}/mapfile
    cd ..
done < genome_mapfile