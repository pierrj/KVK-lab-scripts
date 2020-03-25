#!/bin/bash
#ls *.fastq | awk '{print substr($0,0,5)}' | sort | uniq > mapfile
while read sample;
do
    sbatch --job-name=$sample.spike_in_readcount --export=sample=$sample /global/home/users/pierrj/git/slurm/spikein_readsmapped_loop.slurm
done < mapfile