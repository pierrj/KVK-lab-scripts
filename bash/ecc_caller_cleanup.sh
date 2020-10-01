#!/bin/bash


## clean up script
while read sample
do
cd ${sample}
echo ${sample}
rm ${sample}.mergedandpe.bwamem.bam
rm parallel.confirmed*
rm tmp.*
rm ecccaller_output.${sample}.details.tsv*
rm ecccaller_output.${sample}.bed*
rm parallel.confirmed
rm parallel.plusone.confirmed
rm renamed.filtered.sorted.${sample}.bam
rm renamed.filtered.sorted.${sample}.bam.bai
rm merged.confirmed*
rm lengthfiltered.merged.splitreads.${sample}.renamed.*.bed
rm *.db

# for i in $(seq 0.1 0.1 1.0); do
#     cd $i
#     echo ${i}
#     rm ${sample}_${i}.mergedandpe.bwamem.bam
#     rm parallel.confirmed*
#     rm tmp.*
#     rm ${sample}_${i}.mergedandpe.bwamem.bam
#     rm parallel.confirmed*
#     rm tmp.*
#     rm ecccaller_output.${sample}_${i}.details.tsv*
#     rm ecccaller_output.${sample}_${i}.bed*
#     rm parallel.confirmed
#     rm parallel.plusone.confirmed
#     rm renamed.filtered.sorted.${sample}_${i}.bam
#     rm renamed.filtered.sorted.${sample}_${i}.bam.bai
#     rm *tmp.rarefaction*
#     rm merged.confirmed*
#     rm lengthfiltered.merged.splitreads.${sample}_${i}.renamed.*.bed
#     cd ..
# done

cd ..
done < mapfile




while read sample
do
cd ${sample}
echo ${sample}
rm *.bam
rm *.sam
cd ..
done < mapfile_HP