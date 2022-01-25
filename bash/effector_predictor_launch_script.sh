#!/bin/bash

PROTEOMES_PATH=/global/scratch/users/pierrj/moryzae_pav/gladieux_et_al_2020_data/gladieux_et_al_2021_annotations
MAPFILE=gladieux_et_al_2021_assemblies_mapfile
N_NODES=20


# cd /global/scratch/users/pierrj/moryzae_pav/effector_annotations


if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read genome; do
    echo /global/home/users/pierrj/git/bash/effector_predictor.sh -i ${PROTEOMES_PATH}/${genome}/${genome}_protein.fasta -o ${genome} >> jobqueue
done < $MAPFILE

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

module load seqtk

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.effector_pred --export=node=$node /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done

if [ -f "all_effector_names" ]; then
    rm all_effector_names
fi

while read genome; do
    cat ${genome}.predicted_effectors >> all_effector_names
done < $MAPFILE