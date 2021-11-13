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

## USAGE ##
# this script counts encompassing split reads per gene, averages the counts across samples, then calls eccdna associated and eccdna absent genes
# -n output name or sample name
# -b bedfile of gene locations
# -i file containing list of samples them
# -e file containing list of files containing junction split read locations for each sample

if [ -f "${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    # encompassing split reads are junction split reads that fully encompass a gene
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadspergene
    ## normalize to total number of split reads in the sample
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadspergene > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${ecc_basename}.normalized.splitreadspergene >> ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}

# normalize and average across tech reps, then bio reps
if [ -f "${OUTPUT_NAME}.normalize_table_column" ]; then
    rm ${OUTPUT_NAME}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${OUTPUT_NAME}.normalize_table_column ; done
paste ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn ${OUTPUT_NAME}.normalize_table_column ${SAMPLE_MAPFILE} > ${OUTPUT_NAME}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${OUTPUT_NAME}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${OUTPUT_NAME}.normalized_binned ${OUTPUT_NAME}.normalized.splitreadspergene

# get top ten % of genes for encompassing split reads and call those common genes or eccdna associated genes
top_ten_percent=$(awk '$2!=0' ${OUTPUT_NAME}.normalized.splitreadspergene | wc -l | awk '{print int($1/10)}')
sort -k2,2nr ${OUTPUT_NAME}.normalized.splitreadspergene | head -${top_ten_percent} | awk '{print $1}' > ${OUTPUT_NAME}.common.genes

# genes never found on eccdnas
awk '{ if ($2==0) {print $1}}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${OUTPUT_NAME}.neverfound.genes

awk '{print $1}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${OUTPUT_NAME}.allgenenames