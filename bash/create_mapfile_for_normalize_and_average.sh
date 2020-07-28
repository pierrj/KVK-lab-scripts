#!/bin/bash
while getopts m:t:n: option ### THIS NEEDS TO BE FIXED ####
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
t) FILE=${OPTARG};;
n) NORMALIZE_FILE=${OPTARG};;
esac
done

if [ -f "mapfile_for_normalize_and_average" ]; then
    rm mapfile_for_normalize_and_average
fi

## GET SAMPLE ##

sample=$(head -1 ${MAPFILE} | cut -f1)

start_string_file=$(echo ${FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub(FS $2,"")}1')

end_string_file=$(echo ${FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub($1 FS,"")}1')

start_string_normalize_file=$(echo ${NORMALIZE_FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub(FS $2,"")}1')

end_string_normalize_file=$(echo ${NORMALIZE_FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub($1 FS,"")}1')

while read line; 
do
    sample=$(echo "$line" | cut -f1)
    bio_rep=$(echo "$line" | cut -f2)
    treatment=$(echo "$line" | cut -f3)
    cd ${sample}
    file_path=$(realpath ${start_string_file}${sample}${end_string_file})
    normalize_path=$(realpath ${NORMALIZE_EXTENSION}${sample}.bam)
    cd ..
    echo -e ${file_path}'\t' ${normalize_path}'\t'${sample}'\t'${bio_rep}'\t'${treatment} >> mapfile_for_normalize_and_average
done < ${MAPFILE}