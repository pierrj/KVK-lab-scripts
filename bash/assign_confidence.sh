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

python /global/home/users/pierrj/git/python/merge_eccs.py ${SAMPLE} ${chrom_count}

python /global/home/users/pierrj/git/python/make_coverage_db.py ${SAMPLE}.genomecoverage.filtered.renamed.bed ${chrom_count}

split --number=l/${THREADS} --numeric-suffixes=1 merged.confirmed merged.confirmed

parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_db.py ${SAMPLE} {} ::: $(seq -w 1 ${THREADS})

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed