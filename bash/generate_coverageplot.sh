#!/bin/bash
while getopts m:s:c: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
c) COVFILE=${OPTARG};;
esac
done

awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.normalized_binned.filtered

chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.normalized_binned.filtered > ${SAMPLE}.normalized_binned.filtered.renamed
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT' > tmp.first_row
cat tmp.first_row ${SAMPLE}.normalized_binned.filtered.renamed > ${SAMPLE}.normalized_binned.filtered.renamed.table

y_max=$(awk '{print $3}' ${SAMPLE}.normalized_binned.filtered | sort -nr | head -1 | awk '{print $0/1000}' | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}' | awk '{print $0*1000}')
Rscript --vanilla /global/home/users/pierrj/git/R/make_coverageplot.R ${SAMPLE}.normalized_binned.filtered.renamed.table ${SAMPLE}.manhattanplot.tiff ${y_max}

rm tmp*