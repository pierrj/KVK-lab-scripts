#!/bin/bash

pathfile="/global/scratch/users/pierrj/eccDNA/pipeline_tests/thorsten_script_test" ##add the path to you bam file here
bamfile="guy11.sorted.mergedandpe.G3_1A_bwamem.bam" ##add the name of your bam file here

cd ${pathfile}

#extract split reads mapping in the same orientation
samtools view -f 81 -F 4 ${bamfile} > reverseread1.${bamfile}.sam
samtools view -f 145 -F 4 ${bamfile} > reverseread2.${bamfile}.sam
samtools view -f 65 -F 20 ${bamfile} > forwardread1.${bamfile}.sam
samtools view -f 129 -F 20 ${bamfile} > forwardread2.${bamfile}.sam

#filter reads so that they only map exactly twice and to the same chromosome
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' reverseread1.${bamfile}.sam reverseread1.${bamfile}.sam > samechromosome.exactlytwice.reverseread1.${bamfile}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' reverseread2.${bamfile}.sam reverseread2.${bamfile}.sam > samechromosome.exactlytwice.reverseread2.${bamfile}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' forwardread1.${bamfile}.sam forwardread1.${bamfile}.sam > samechromosome.exactlytwice.forwardread1.${bamfile}.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' forwardread2.${bamfile}.sam forwardread2.${bamfile}.sam > samechromosome.exactlytwice.forwardread2.${bamfile}.sam

#combine all into a single file and convert back to bam file
cat samechromosome.exactlytwice.reverseread1.${bamfile}.sam samechromosome.exactlytwice.reverseread2.${bamfile}.sam samechromosome.exactlytwice.forwardread1.${bamfile}.sam samechromosome.exactlytwice.forwardread2.${bamfile}.sam > samechromosome.exactlytwice.all.${bamfile}.sam
samtools view -b -h <(cat <(samtools view -H ${bamfile}) samechromosome.exactlytwice.all.${bamfile}.sam) > splitreads.${bamfile}