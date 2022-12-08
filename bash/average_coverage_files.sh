#!/bin/bash
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
    bedtools coverage -sorted -counts -a $BED_FILE -g ${GENOME_CHROMSIZES} -b $density_file | awk '{print $(NF)}' > ${basename_density_file}.${LABEL}.coverage
    paste ${basename_density_file}.${LABEL}.coverage ${LABEL}.lengths | awk '{print $1/($2)}' > ${basename_density_file}.${LABEL}.coverage.RPK
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${basename_density_file}.${LABEL}.coverage.RPK)
    awk -v N=$sum_rpks '{print $1/N}' ${basename_density_file}.${LABEL}.coverage.RPK > ${basename_density_file}.${LABEL}.${OUTPUT_NAME}.coverage.normalized
done < $BAMFILE_MAPFILE
paste *.${LABEL}.${OUTPUT_NAME}.coverage.normalized | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${LABEL}.${OUTPUT_NAME}.txt