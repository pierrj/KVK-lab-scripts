#!/bin/bash

PROTEOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/all_proteomes_corrected/
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/proteomes_mapfile
N_NODES=20


cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/predicted_effectors

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read proteome; do
    echo /global/home/users/pierrj/git/bash/effector_predictor.sh -i ${PROTEOMES_PATH}/${proteome} -o ${proteome}.effectors_pred >> jobqueue
done < ${MAPFILE}

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.effector_pred --export=ALL,node=$node /global/home/users/pierrj/git/slurm/gnu_parallel_multinode_v2.slurm
done

if [ -f "all_effector_names" ]; then
    rm all_effector_names
fi

while read genome; do
    cat ${genome}.effectors_pred.predicted_effectors >> all_effector_names
done < $MAPFILE