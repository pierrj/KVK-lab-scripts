#!/bin/bash
ls *.fastq | awk '{print substr($0,0,5)}' | sort | uniq > mapfile
while read sample;
do
    sbatch --job-name=$sample.ecc_caller --export=sample=$sample /global/home/users/pierrj/git/slurm/ecc_callerv3_pureculture.slurm
done < mapfile