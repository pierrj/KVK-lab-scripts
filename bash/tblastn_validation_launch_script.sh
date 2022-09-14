#!/bin/bash

PROJECT_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all
DATE=Sep08

cd $PROJECT_DIR

LOST_GENOME_DIR=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate # where the assemblies are
LOST_OG_DIR=${PROJECT_DIR}/orthofinder/Results_${DATE}/WorkingDirectory/OrthoFinder/Results_${DATE}/Orthogroup_Sequences # where the OG_protein fastas are
E_VALUE=1e-10
PIDENT=55
QUERY_COV=55
HIT_COUNT=2
N_NODES=40
OUTPUT_FILE=${PROJECT_DIR}/pav_table
BLAST_DB=${PROJECT_DIR}/all_ogs_seqs.fasta
ABSENCES_FILE=${PROJECT_DIR}/absences_to_validate.tsv ## location of table of absences to validate

conda activate /global/scratch/users/pierrj/conda_envs/orthofinder

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/make_single_file_from_og_dir.py ${LOST_OG_DIR} ${BLAST_DB}

makeblastdb -in ${BLAST_DB} -dbtype prot

tr -d '\015' <${ABSENCES_FILE} > ${ABSENCES_FILE}.fixed

if [ -d "${PROJECT_DIR}/pav_validation" ]; then
    rm -r ${PROJECT_DIR}/pav_validation
fi

mkdir ${PROJECT_DIR}/pav_validation

cd ${PROJECT_DIR}/pav_validation

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}.fa -g ${LOST_GENOME_DIR}/${LOST_GENOME} \
        -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -d ${BLAST_DB} >> jobqueue
done < ${ABSENCES_FILE}.fixed

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