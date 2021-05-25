#!/bin/bash
while getopts m:s:t:b:r: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};;
r) CONFIRMED_SPLITREADS=${OPTARG};;
esac
done



chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
samtools view -H ${FILTERED_BAMFILE} | awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{ if ($2 ~ /^SN/ && substr($2, 4) in a) {print $1, "SN:"a[substr($2,4)], $3} else {print $0}}' tmp.chrom_count_and_names - |\
    samtools reheader - ${FILTERED_BAMFILE} > renamed.filtered.sorted.${SAMPLE}.bam

samtools index renamed.filtered.sorted.${SAMPLE}.bam

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${CONFIRMED_SPLITREADS} > parallel.plusone.confirmed
awk -v OFS='\t' '{print $1-1, $2, $3}' parallel.plusone.confirmed > parallel.confirmed

## this is a little messy

sed 's/[[:space:]]*$//' parallel.confirmed | sort -k1,1n -k2,2n  | uniq -c | awk -v OFS='\t' '{print $2, $3, $4, $1}' > merged.confirmed

shuf merged.confirmed > shuf.merged.confirmed

split --number=l/${THREADS} --numeric-suffixes=1 shuf.merged.confirmed merged.confirmed

parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_nodb_variablesrs.py ${SAMPLE} {} renamed.filtered.sorted.${SAMPLE}.bam 2 4 ::: $(seq -w 1 ${THREADS})

# put parallel chunks back together
cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv
cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed
cat $(find . -maxdepth 1 -name "ecccaller_output.splitreads.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.splitreads.${SAMPLE}.bed

# rename output files to original chrom/scaffold names
paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.details.tsv
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.splitreads.${SAMPLE}.bed > ecccaller_output.splitreads.${SAMPLE}.renamed.bed

sort -k1,1 -k2,2n ecccaller_output.${SAMPLE}.renamed.details.tsv > ${SAMPLE}.ecc_caller_out.details.txt
sort -k1,1 -k2,2n ecccaller_output.${SAMPLE}.renamed.bed > ${SAMPLE}.ecc_caller_out.genomebrowser.bed
sort -k1,1 -k2,2n ecccaller_output.splitreads.${SAMPLE}.renamed.bed > ${SAMPLE}.ecc_caller_out.splitreads.bed

# clean up tmp files
rm ecccaller_output.${SAMPLE}.details.tsv*
rm ecccaller_output.${SAMPLE}.bed*
rm ecccaller_output.splitreads.${SAMPLE}.bed*
rm parallel.confirmed
rm parallel.plusone.confirmed
rm renamed.filtered.sorted.${SAMPLE}.bam
rm renamed.filtered.sorted.${SAMPLE}.bam.bai
rm merged.confirmed*
rm shuf.merged.confirmed
rm tmp.*
rm ecccaller_output.${SAMPLE}.renamed.details.tsv
rm ecccaller_output.${SAMPLE}.renamed.bed
rm ecccaller_output.splitreads.${SAMPLE}.renamed.bed