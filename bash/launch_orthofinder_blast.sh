#!/bin/bash

N_NODES=50

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast

conda activate /global/scratch/users/pierrj/conda_envs/orthofinder

orthofinder -op -S diamond_ultra_sens -f all_proteomes_corrected -o orthofinder | grep "diamond blastp" > jobqueue

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split -a 3 --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%03g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio2 --ntasks-per-node 24 --job-name=$node.blast --export=node=$node /global/home/users/pierrj/git/slurm/orthofinder_blast.slurm
done