#!/bin/bash
while getopts m:s:t:c: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
esac
done

# while read sample; 
# do
#     cd ${sample}
#     normalize_factor=$(samtools view -c -F 4 -F 2048 sorted.MG8.mergedandpe.${sample}_bwamem.bam | awk '{print $1/1000000}' )
#     awk -v N=$normalize_factor '{sum+=$3} NR%100==0 {print sum/100/N; sum =0}' genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/${sample}.${table}
#     cd ..
# done < mapfile