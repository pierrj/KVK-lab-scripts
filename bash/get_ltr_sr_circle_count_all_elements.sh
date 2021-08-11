#!/bin/bash
while getopts b:m:s:o: option
do
case "${option}"
in
b) BEDFILE=${OPTARG};;
m) MAPFILE=${OPTARG};;
s) SPLIT_READ_FILE=${OPTARG};;
o) OUTPUTNAME=${OPTARG};;
esac
done

if [ -f "${OUTPUTNAME}.sr_count_per_element" ]; then
    rm ${OUTPUTNAME}.sr_count_per_element
fi

if [ -f "${OUTPUTNAME}.ltr_sr_cov_perfeature" ]; then
    rm ${OUTPUTNAME}.ltr_sr_cov_perfeature
fi

if [ -f "${OUTPUTNAME}.read_cov_perfeature" ]; then
    rm ${OUTPUTNAME}.read_cov_perfeature
fi


## add to file for all elements

while read element
do
/global/home/users/pierrj/git/bash/get_ltr_sr_circle_count_per_element.sh -b ${BEDFILE} -e ${element} -s ${SPLIT_READ_FILE} -o ${OUTPUTNAME} >> ${OUTPUTNAME}.sr_count_per_element
cat ${ELEMENT}.${OUTPUTNAME}.ltr_sr_cov_perfeature >> ${OUTPUTNAME}.ltr_sr_cov_perfeature
cat ${ELEMENT}.${OUTPUTNAME}.read_cov_perfeature >> ${OUTPUTNAME}.read_cov_perfeature
done < ${MAPFILE}