#!/bin/bash
while getopts s:f:l:e: option
do
case "${option}"
in
s) SAMPLE=${OPTARG};;
f) SPLITREAD_FILE=${OPTARG};;
l) LTR_TE_FILE=${OPTARG};;
e) ELEMENT_MAPFILE=${OPTARG};;
esac
done

bedtools intersect -wao -a ${SPLITREAD_FILE} -b ${LTR_TE_FILE} | \
     awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
     awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > ${SAMPLE}.ltr_eccs

ltr_split_read_total=$(wc -l ${SAMPLE}.ltr_eccs | awk '{print $1}')

if [ -f "${OUTPUTNAME}.sr_count_per_element" ]; then
    rm ${OUTPUTNAME}.sr_count_per_element
fi

while read ELEMENT
do
grep ${ELEMENT} ${LTR_TE_FILE} > ${ELEMENT}.loc.gff
bedtools intersect -wao -a ${SPLITREAD_FILE} -b ${ELEMENT}.loc.gff | \
     awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
     awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > ${SAMPLE}.${ELEMENT}.ltr_eccs
sr_count=$(wc -l ${SAMPLE}.${ELEMENT}.ltr_eccs | awk '{print $1}')
echo -e ${ELEMENT}'\t'${sr_count} >> ${OUTPUTNAME}.sr_count_per_element
done < ${ELEMENT_MAPFILE}