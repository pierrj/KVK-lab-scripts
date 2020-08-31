#!/bin/bash
while getopts b:m:s:o:f: option
do
case "${option}"
in
b) BEDFILE=${OPTARG};;
m) MAPFILE=${OPTARG};;
s) SPLIT_READ_FILE=${OPTARG};;
o) OUTPUTNAME=${OPTARG};;
f) FILTERED_BAMFILE=${OPTARG};;
esac
done

if [ -f "${OUTPUTNAME}.sr_count_per_element" ]; then
    rm ${OUTPUTNAME}.sr_count_per_element
fi

while read element
do
/global/home/users/pierrj/git/bash/get_ltr_sr_circle_count_per_element.sh -b ${BEDFILE} -e ${element} -s ${SPLIT_READ_FILE} -o ${OUTPUTNAME} -f ${FILTERED_BAMFILE} >> ${OUTPUTNAME}.sr_count_per_element
done < ${MAPFILE}