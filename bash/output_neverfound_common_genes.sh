#!/bin/bash
while getopts n:b:i:e: option
do
case "${option}"
in
n) OUTPUT_NAME=${OPTARG};;
b) GENE_BEDFILE=${OPTARG};;
i) SAMPLE_MAPFILE=${OPTARG};;
e) ECCDNA_MAPFILE=${OPTARG};;
esac
done

if [ -f "${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadspergene > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${ecc_basename}.normalized.splitreadspergene >> ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}

if [ -f "${OUTPUT_NAME}.normalize_table_column" ]; then
    rm ${OUTPUT_NAME}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${OUTPUT_NAME}.normalize_table_column ; done
paste ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn ${OUTPUT_NAME}.normalize_table_column ${SAMPLE_MAPFILE} > ${OUTPUT_NAME}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${OUTPUT_NAME}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${OUTPUT_NAME}.normalized_binned ${OUTPUT_NAME}.normalized.splitreadspergene

## maybe something besides 1000 here?
# maybe top 10% of found genes like this
top_ten_percent=$(awk '$2!=0' ${OUTPUT_NAME}.normalized.splitreadspergene | wc -l | awk '{print int($1/10)}')
sort -k2,2nr ${OUTPUT_NAME}.normalized.splitreadspergene | head -${top_ten_percent} | awk '{print $1}' > ${OUTPUT_NAME}.common.genes

awk '{ if ($2==0) {print $1}}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${OUTPUT_NAME}.neverfound.genes

awk '{print $1}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${OUTPUT_NAME}.allgenenames