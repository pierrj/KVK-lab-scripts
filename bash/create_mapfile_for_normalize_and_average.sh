#!/bin/bash
while getopts m:t:n: option ### THIS NEEDS TO BE FIXED ####
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
t) FILE_EXTENSION=${OPTARG};;
n) NORMALIZE_EXTENSION=${OPTARG};;
esac
done

if [ -f "mapfile_for_normalize_and_average" ]; then
    rm mapfile_for_normalize_and_average
fi

filtered.sorted.RC_1A.bam

while read line; 
do
    sample=$(echo "$line" | cut -f1)
    bio_rep=$(echo "$line" | cut -f2)
    treatment=$(echo "$line" | cut -f3)
    cd ${sample}
    file_path=$(realpath ${sample}${FILE_EXTENSION})
    normalize_path=$(realpath ${NORMALIZE_EXTENSION}${sample}.bam)
    cd ..
    echo -e ${file_path}'\t' ${normalize_path}'\t'${bio_rep}'\t'${treatment} >> mapfile_for_normalize_and_average
done < ${MAPFILE}