#!/bin/bash
while getopts m:s:t:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};;
esac
done



chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
samtools view -H ${FILTERED_BAMFILE} | awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{ if ($2 ~ /^SN/ && substr($2, 4) in a) {print $1, "SN:"a[substr($2,4)], $3} else {print $0}}' tmp.chrom_count_and_names - |\
    samtools reheader - ${FILTERED_BAMFILE} > renamed.filtered.sorted.${SAMPLE}.bam

samtools index renamed.filtered.sorted.${SAMPLE}.bam

python /global/home/users/pierrj/git/python/merge_eccs.py ${SAMPLE} ${chrom_count}

split --number=l/${THREADS} --numeric-suffixes=1 merged.confirmed merged.confirmed

parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_nodb.py ${SAMPLE} {} renamed.filtered.sorted.${SAMPLE}.bam ::: $(seq -w 1 ${THREADS})

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed

paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed