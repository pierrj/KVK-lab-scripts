#!/bin/bash

cd /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/og_unannotated_full_run

LOST_GENOME_DIR=/global/scratch/users/pierrj/fungap_runs/gladieux_all/genomes_to_annotate # where the assemblies are
LOST_OG_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022/orthofinder_1_31_2022/Results_Feb03/WorkingDirectory/OrthoFinder/Results_Feb03/Orthogroup_Sequences # where the OG_protein fastas are
E_VALUE=1e-10
PIDENT=55
QUERY_COV=55
HIT_COUNT=2
N_NODES=20
OUTPUT_FILE=tblastn_pav_validation_full_run_2_7_2022
GENE_GFF=/global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/CH0043_merge.gff3
GENE_OVERLAP=90


if [ -f "jobqueue" ]; then
    rm jobqueue
fi


while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}.fa -g ${LOST_GENOME_DIR}/${LOST_GENOME} \
        -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -f ${GENE_GFF} -b ${GENE_OVERLAP} >> jobqueue
done < /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/absences_to_validate

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

if [ -f "${OUTPUT_FILE}" ]; then
    rm ${OUTPUT_FILE}
fi

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    echo $node
    sbatch --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE /global/home/users/pierrj/git/slurm/gnu_parallel_multinode.slurm
done