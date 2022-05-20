#!/bin/bash

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pav_validation

LOST_GENOME_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/assemblies # where the assemblies are
LOST_OG_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/orthofinder_4_27_2022/Results_Apr28/WorkingDirectory/OrthoFinder/Results_Apr28/Orthogroup_Sequences # where the OG_protein fastas are
E_VALUE=1e-10
PIDENT=55
QUERY_COV=55
HIT_COUNT=2
N_NODES=40
OUTPUT_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pav_table
BLAST_DB=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/all_ogs_seqs.fasta
ABSENCES_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/absences_to_validate.txt ## location of table of absences to validate


if [ -f "jobqueue" ]; then
    rm jobqueue
fi


while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}.fa -g ${LOST_GENOME_DIR}/${LOST_GENOME} \
        -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -d ${BLAST_DB} >> jobqueue
done < ${ABSENCES_FILE}

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

if [ -f "${OUTPUT_FILE}" ]; then
    rm ${OUTPUT_FILE}
fi

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio2 --ntasks-per-node 24 --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done

## after everything is done

awk -v OFS='\t' '{ if ($3 == "yes") {print $1, $2, $3} else {print $1, $2, "no"}}' ${OUTPUT_FILE} > ${OUTPUT_FILE}.simplified



## ugh test old version...

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pav_validation_old_version

LOST_GENOME_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/assemblies # where the assemblies are
LOST_OG_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/orthofinder_4_27_2022/Results_Apr28/WorkingDirectory/OrthoFinder/Results_Apr28/Orthogroup_Sequences # where the OG_protein fastas are
E_VALUE=1e-10
PIDENT=55
QUERY_COV=55
HIT_COUNT=2
N_NODES=40
OUTPUT_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pav_table_old_version
ABSENCES_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/absences_to_validate.txt ## location of table of absences to validate

if [ -f "jobqueue" ]; then
    rm jobqueue
fi


while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation_dec_27_version.sh -l ${LOST_OG_DIR}/${LOST_OG}.fa -g ${LOST_GENOME_DIR}/${LOST_GENOME} \
        -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} >> jobqueue
done < ${ABSENCES_FILE}

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

if [ -f "${OUTPUT_FILE}" ]; then
    rm ${OUTPUT_FILE}
fi

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio2 --ntasks-per-node 24 --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done