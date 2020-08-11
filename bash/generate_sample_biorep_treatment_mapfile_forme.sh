#!/bin/bash
while getopts m: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
esac
done

if [ -f "sample_mapfile" ]; then
    rm sample_mapfile
fi

while read sample; 
do
    bio_rep=$(echo ${sample} | cut -c-4)
    treatment=$(echo ${sample} | cut -c-2)
    echo -e ${sample}'\t'${bio_rep}'\t'${treatment} >> sample_mapfile
done < ${MAPFILE}