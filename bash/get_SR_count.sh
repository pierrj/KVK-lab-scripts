#!/bin/bash
samtools view -f 81 -F 4 sorted.mergedandpe.${1}_bwamem.bam ${2} > tmp.reverseread1.mergedandpe.${2}_bwamem.sam
samtools view -f 145 -F 4 sorted.mergedandpe.${1}_bwamem.bam ${2} > tmp.reverseread2.mergedandpe.${2}_bwamem.sam
samtools view -f 65 -F 20 sorted.mergedandpe.${1}_bwamem.bam ${2} > tmp.forwardread1.mergedandpe.${2}_bwamem.sam
samtools view -f 129 -F 20 sorted.mergedandpe.${1}_bwamem.bam ${2} > tmp.forwardread2.mergedandpe.${2}_bwamem.sam

awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${2}_bwamem.sam tmp.reverseread1.mergedandpe.${2}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${2}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${2}_bwamem.sam tmp.reverseread2.mergedandpe.${2}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${2}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${2}_bwamem.sam tmp.forwardread1.mergedandpe.${2}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${2}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${2}_bwamem.sam tmp.forwardread2.mergedandpe.${2}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${2}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${2}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${2}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${2}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${2}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam


awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam > qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam > exactlytwice.qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam
samtools view -b -h <(cat <(samtools view -H mergedandpe.${1}_bwamem.bam) exactlytwice.qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${2}_bwamem.bam) > SRs.mergedandpe.${2}_bwamem.bam
bedtools bamtobed -i SRs.mergedandpe.${2}_bwamem.bam | sort -k 4,4 -k 2,2 > SRs.mergedandpe.${2}_bwamem.bed

rm tmp.*
wc -l SRs.mergedandpe.${2}_bwamem.bed | awk '{print $1/2}'