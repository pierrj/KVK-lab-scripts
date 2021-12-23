#!/bin/bash

LOST_GENOME_DIR=
LOST_OG_DIR=
E_VALUE=
PIDENT=
QUERY_COV=
HIT_COUNT=
N_NODES=20
OUTPUT_FILE=tblastn_validation_out


if [ -f "jobqueue" ]; then
    rm jobqueue
fi


while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}_protein.fasta -g ${LOST_GENOME_DIR}/${LOST_GENOME} -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} >> jobqueue
done < absences_to_validate.tsv

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