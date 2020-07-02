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

ipcluster start -n ${THREADS} --cluster-id="cluster-id-${SAMPLE}" --profile=pierrj &
sleep 300
ipython /global/home/users/pierrj/git/python/ecc_caller_anygenome_assignconfidence_slow.py ${SAMPLE}.genomecoverage.filtered.renamed.bed ${SAMPLE} ${chrom_count} ${SORTED_BAMFILE}
ipcluster stop --cluster-id="cluster-id-${SAMPLE}"

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed