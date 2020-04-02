#!/bin/bash
#ls *.fastq | awk '{print substr($0,0,5)}' | sort | uniq > mapfile
while read sample;
if [ ${sample:0:2} = "IF" ]
do
    sbatch --job-name=$sample.unmapped_reads --export=sample=$sample /global/home/users/pierrj/git/slurm/blast_unmapped_reads_loop.slurm
done < mapfile