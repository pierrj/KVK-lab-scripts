#!/bin/bash
while getopts p:n:c:s:g: option
do
case "${option}"
in
p) PANNZER_INPUT=${OPTARG};;
n) OUTPUT_NAME=${OPTARG};;
c) SCORE_CUTOFF=${OPTARG};;
s) SIGNIFICANT_GENES=${OPTARG};;
g) GO_TERM=${OPTARG};;
esac
done

python /global/home/users/pierrj/git/python/parse_pannzer_output.py ${PANNZER_INPUT} ${SCORE_CUTOFF} ${OUTPUT_NAME}.filtered.GO.txt

awk -v OFS='\t' '{print $0}' ${OUTPUT_NAME}.filtered.GO.txt > ${OUTPUT_NAME}.filtered.GO.awk.txt

Rscript --vanilla /global/home/users/pierrj/git/R/go_enrichment.R ${OUTPUT_NAME}.filtered.GO.awk.txt ${SIGNIFICANT_GENES} ${GO_TERM} ${OUTPUT_NAME}.results.txt

rm ${OUTPUT_NAME}.filtered.GO.awk.txt

rm ${OUTPUT_NAME}.filtered.GO.txt