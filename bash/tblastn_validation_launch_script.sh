#!/bin/bash

LOST_GENOME_DIR=/global/scratch/users/pierrj/moryzae_pav/gladieux_et_al_2020_data/gladieux_et_al_2021_assemblies
LOST_OG_DIR=/global/scratch/users/pierrj/moryzae_pav/orthogroup_cons/orthogroup_protein_out
E_VALUE=1e-10
PIDENT=55
QUERY_COV=55
HIT_COUNT=2
N_NODES=20
OUTPUT_FILE=tblastn_pav_validation_gladieux_1_3_2022


if [ -f "jobqueue" ]; then
    rm jobqueue
fi


while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}_protein.fasta -g ${LOST_GENOME_DIR}/${LOST_GENOME}.fasta -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} >> jobqueue
done < absences_to_validate

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

if [ -f "${OUTPUT_FILE}" ]; then
    rm ${OUTPUT_FILE}
fi

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done