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

samtools view -b -q 1 ${SORTED_BAMFILE} $(cat ${MAPFILE} | tr "\n" " ") > filtered.${SORTED_BAMFILE}

##MAKE SURE THIS GETS FIXED CURRENTLY THIS IS JUST FOR TESTING
samtools view -f 81 -F 4 filtered.${SORTED_BAMFILE}> tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 filtered.${SORTED_BAMFILE} > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 filtered.${SORTED_BAMFILE} > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 filtered.${SORTED_BAMFILE} > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
## samtools view -f 16 -F 5 filtered.${SORTED_BAMFILE} > tmp.reversemerged.mergedandpe.${sample}_bwamem.sam
## samtools view -F 21 filtered.${SORTED_BAMFILE} > tmp.forwardmerged.mergedandpe.${sample}_bwamem.sam

awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

samtools view guy11.sorted.mergedandpe.${sample}_bwamem.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) - | samtools view -b -h - > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed
mv tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed > outwardfacing.mergedandpe.${sample}_bwamem.bed
