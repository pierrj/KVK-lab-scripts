#!/bin/bash
while getopts m:s:t:c: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
c) COVFILE=${OPTARG};;
esac
done

awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.genomecoverage.filtered.bed
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.genomecoverage.filtered.bed > ${SAMPLE}.genomecoverage.filtered.renamed.bed

normalize_factor=$(samtools view -c -F 4 -F 2048 ${SAMPLE}.sorted.mergedandpe.bwamem.bam | awk '{print $1/1000000}' )
awk -v N=$normalize_factor -v OFS='\t' '{sum+=$3} NR%100==0 {print $1, $2, sum/100/N; sum =0}' ${SAMPLE}.genomecoverage.bed > ${SAMPLE}.genomecoverage.normalized.bed
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT' > tmp.first_row
cat tmp.first_row ${SAMPLE}.genomecoverage.normalized.bed > ${SAMPLE}.genomecoverage.table

y_max=$(awk '{print $3}' ${SAMPLE}.genomecoverage.normalized.bed | sort -nr | head -1 | awk '{print $0/1000}' | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}' | awk '{print $0*1000}')
Rscript --vanilla /global/home/users/pierrj/git/R/make_coverageplot.R ${SAMPLE}.genomecoverage.table ${SAMPLE}.genomecoverage.tiff ${y_max}

rm tmp*