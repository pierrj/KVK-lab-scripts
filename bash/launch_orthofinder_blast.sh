#!/bin/bash

N_NODES=100

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes

conda activate /global/scratch/users/pierrj/conda_envs/orthofinder

orthofinder -op -S diamond_ultra_sens -f all_proteomes_corrected -o orthofinder_4_27_2022 | tail -n +184 > jobqueue

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split -a 3 --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%03g" 2 ${N_NODES})
do
    sbatch --job-name=$node.blast --export=node=$node /global/home/users/pierrj/git/slurm/orthofinder_blast.slurm
done

for node in $(seq -f "%03g" 1 ${N_NODES})
do
    echo $node
done