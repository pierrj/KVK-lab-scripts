#!/bin/bash
samtools view -f 81 -F 4 mergedandpe.${1}_bwamem.bam > tmp.reverseread1.mergedandpe.${1}_bwamem.sam
samtools view -f 145 -F 4 mergedandpe.${1}_bwamem.bam > tmp.reverseread2.mergedandpe.${1}_bwamem.sam
samtools view -f 65 -F 20 mergedandpe.${1}_bwamem.bam > tmp.forwardread1.mergedandpe.${1}_bwamem.sam
samtools view -f 129 -F 20 mergedandpe.${1}_bwamem.bam > tmp.forwardread2.mergedandpe.${1}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${1}_bwamem.sam tmp.reverseread1.mergedandpe.${1}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${1}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${1}_bwamem.sam tmp.reverseread2.mergedandpe.${1}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${1}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${1}_bwamem.sam tmp.forwardread1.mergedandpe.${1}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${1}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${1}_bwamem.sam tmp.forwardread2.mergedandpe.${1}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${1}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${1}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${1}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${1}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${1}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${1}_bwamem.bam

samtools view mergedandpe.${1}_bwamem.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H mergedandpe.${1}_bwamem.bam) - | samtools view -b -h - > tmp.outwardfacing.mergedandpe.${1}_bwamem.bam
bedtools bamtobed -i tmp.outwardfacing.mergedandpe.${1}_bwamem.bam | sort -k 4,4 > outwardfacing.mergedandpe.${1}_bwamem.bed

samtools sort mergedandpe.${1}_bwamem.bam -o tmp.sorted.mergedandpe.${1}_bwamem.bam
bedtools genomecov -d -ibam tmp.sorted.mergedandpe.${1}_bwamem.bam > genomecoverage.mergedandpe.${1}_bwamem.bed

rm tmp.*