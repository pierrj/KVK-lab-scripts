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
ABSENCES_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/absences_to_validate.tsv ## location of table of absences to validate


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
    sbatch -p savio --ntasks-per-node 20 --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done






cd /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/test_blast_id

LOST_GENOME=/global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/CH0043.fasta
E_VALUE=1e-10
HIT_COUNT=2
QUERY_COV=55
PIDENT=55
BLAST_DB=/global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/all_proteins_with_og_names_no_CH0043.fasta
N_NODES=40
OUTPUT_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/test_blast_id/optimization_out_5_9_22_wout_ch0043_in_blastdb

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read og; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate/orthogroup_protein_out_no_CH0043/${og}_protein.fasta \
    -g ${LOST_GENOME} -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -d ${BLAST_DB} >> jobqueue
done < /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate/og_list_all_presences_CH0043

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

if [ -f "${OUTPUT_FILE}" ]; then
    rm ${OUTPUT_FILE}
fi

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio --ntasks-per-node 20 --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done