#!/bin/bash
#ls *.fastq | awk '{print substr($0,0,5)}' | sort | uniq > mapfile
while read sample;
do
    sbatch --job-name=$sample.spike_in_SRs --export=sample=$sample /global/home/users/pierrj/git/slurm/spikein_SRs_loop.slurm
done < mapfile