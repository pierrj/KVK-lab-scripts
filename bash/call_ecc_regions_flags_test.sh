#!/bin/bash
while getopts m:s:t:c:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
c) COVFILE=${OPTARG};;
b) SORTED_BAMFILE=${OPTARG};;
esac
done

### THIS VERSION IS A TEST TO MAKE SURE THE FLAG OPTIONS REPLICATE THE NO FLAGS PIPELINE

samtools view -b -q 1 ${SORTED_BAMFILE} $(cat ${MAPFILE} | tr "\n" " ") > filtered.${SAMPLE}.bam

samtools view -f 81 -F 4 filtered.${SAMPLE}.bam > tmp.reverseread1.${SAMPLE}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.${SAMPLE}.sam tmp.reverseread1.${SAMPLE}.sam > tmp.samechromosome.exactlytwice.reverseread1.${SAMPLE}.sam

samtools view -f 145 -F 4 filtered.${SAMPLE}.bam > tmp.reverseread2.${SAMPLE}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.${SAMPLE}.sam tmp.reverseread2.${SAMPLE}.sam > tmp.samechromosome.exactlytwice.reverseread2.${SAMPLE}.sam

samtools view -f 65 -F 20 filtered.${SAMPLE}.bam > tmp.forwardread1.${SAMPLE}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.${SAMPLE}.sam tmp.forwardread1.${SAMPLE}.sam > tmp.samechromosome.exactlytwice.forwardread1.${SAMPLE}.sam

samtools view -f 129 -F 20 filtered.${SAMPLE}.bam > tmp.forwardread2.${SAMPLE}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.${SAMPLE}.sam tmp.forwardread2.${SAMPLE}.sam > tmp.samechromosome.exactlytwice.forwardread2.${SAMPLE}.sam

cat tmp.samechromosome.exactlytwice.reverseread1.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.reverseread2.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.forwardread1.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.forwardread2.${SAMPLE}.sam > samechromosome.exactlytwice.all.mergedandpe.${SAMPLE}_bwamem.sam

samtools view filtered.${SAMPLE}.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H ${SORTED_BAMFILE}) - | samtools view -b -h - > tmp.outwardfacing.${SAMPLE}.bam
bedtools bamtobed -i tmp.outwardfacing.${SAMPLE}.bam | sort -k 4,4 > tmp.outwardfacing.${SAMPLE}.bed
mv tmp.outwardfacing.${SAMPLE}.bed tmp.outwardfacing.${SAMPLE}.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.${SAMPLE}.bed.old > tmp.outwardfacing.${SAMPLE}.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.${SAMPLE}.bed.old.trimmed tmp.outwardfacing.${SAMPLE}.bed.old.trimmed > outwardfacing.${SAMPLE}.bed

awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.genomecoverage.filtered.bed
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.genomecoverage.filtered.bed > ${SAMPLE}.genomecoverage.filtered.renamed.bed
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names outwardfacing.${SAMPLE}.bed > outwardfacing.${SAMPLE}.renamed.bed

ipcluster start -n ${THREADS} --cluster-id="cluster-id-${SAMPLE}" &
sleep 300
ipython /global/home/users/pierrj/git/python/ecc_caller_anygenome.py samechromosome.exactlytwice.all.mergedandpe.${SAMPLE}_bwamem.sam outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE}.genomecoverage.filtered.renamed.bed ${SAMPLE} ${chrom_count} ${SORTED_BAMFILE}
ipcluster stop --cluster-id="cluster-id-${SAMPLE}"

paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed
awk -v OFS='\t' '{print $1+1, $2, $3}' parallel.confirmed > parallel.confirmed.plusone
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count parallel.confirmed.plusone > ${SAMPLE}.confirmedsplitreads.bed

rm tmp.*