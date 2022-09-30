#!/bin/bash

N_NODES=50

PROJECT_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap

cd $PROJECT_DIR

conda activate /global/scratch/users/pierrj/conda_envs/orthofinder

# orthofinder -op -S diamond_ultra_sens -f all_proteomes_corrected -o orthofinder | grep "diamond blastp" > jobqueue

# orthofinder -op -S diamond_ultra_sens -f all_proteomes_corrected_w_70_15 -o orthofinder_w_70_15 | grep "diamond blastp" > jobqueue

orthofinder -op -S diamond_ultra_sens -f all_proteomes_corrected_w_guy11 -o orthofinder_w_guy11 | grep "diamond blastp" > jobqueue

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split -a 3 --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%03g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio2 --ntasks-per-node 24 --job-name=$node.blast --export=ALL,node=$node /global/home/users/pierrj/git/slurm/orthofinder_blast.slurm
done