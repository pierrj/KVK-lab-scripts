#!/bin/bash
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
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
while getopts b:l:m:g:o: option
do
case "${option}"
in
b) BED_FILE=${OPTARG};;
l) LABEL=${OPTARG};;
m) BAMFILE_MAPFILE=${OPTARG};;
g) GENOME_CHROMSIZES=${OPTARG};;
o) OUTPUT_NAME=${OPTARG};;
esac
done

while read density_file; do
    basename_density_file=$(basename $density_file)
    # calculate coverage for each bed entry
    bedtools coverage -sorted -counts -a $BED_FILE -g ${GENOME_CHROMSIZES} -b $density_file | awk '{print $(NF)}' > ${basename_density_file}.${LABEL}.coverage
    # paste bed lengths (calculated before this script), then normalize by length to get RPK values
    paste ${basename_density_file}.${LABEL}.coverage ${LABEL}.lengths | awk '{print $1/($2)}' > ${basename_density_file}.${LABEL}.coverage.RPK
    # sum RPKS
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${basename_density_file}.${LABEL}.coverage.RPK)
    # calculated per million of RPKs for each bed entry
    awk -v N=$sum_rpks '{print $1/N}' ${basename_density_file}.${LABEL}.coverage.RPK > ${basename_density_file}.${LABEL}.${OUTPUT_NAME}.coverage.normalized
done < $BAMFILE_MAPFILE
# paste and average all files together
paste *.${LABEL}.${OUTPUT_NAME}.coverage.normalized | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${LABEL}.${OUTPUT_NAME}.txt