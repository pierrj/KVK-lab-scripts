#!/bin/bash
while getopts m:s:c:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
c) COVFILE=${OPTARG};;
b) BEDFILE=${OPTARG};;
esac
done

awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.normalized_binned.filtered

chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.normalized_binned.filtered > ${SAMPLE}.normalized_binned.filtered.renamed

## for highlighting, probably to be made optional
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${BEDFILE} > renamed.${SAMPLE}.highlight_locs.bed
awk -v OFS='\t' '{print $1, $2, $2+99, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > binned.${SAMPLE}.normalized_binned.filtered.renamed.bed
bedtools intersect -wa -f 0.1 -a binned.${SAMPLE}.normalized_binned.filtered.renamed.bed -b renamed.${SAMPLE}.highlight_locs.bed | cut -f 4 > snps.tohighlight
awk -v OFS='\t' '{print $1, $2, $3, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT''\t''SNP' > tmp.first_row
cat tmp.first_row ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table

y_max=$(awk '{print $3}' ${SAMPLE}.normalized_binned.filtered | sort -nr | head -1 | awk '{print $0/1000}' | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}' | awk '{print $0*1000}')
Rscript --vanilla /global/home/users/pierrj/git/R/make_coverageplot_highlight.R ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table ${SAMPLE}.manhattanplot.tiff ${y_max} snps.tohighlight

rm tmp*