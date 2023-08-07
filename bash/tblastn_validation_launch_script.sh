#!/bin/bash
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

PROJECT_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_anne_subset_re
DATE=out

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

# concatenate all orthogroup sequences together into one file and make blast db
/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/make_single_file_from_og_dir.py ${LOST_OG_DIR} ${BLAST_DB}

makeblastdb -in ${BLAST_DB} -dbtype prot

# fix line endings
tr -d '\015' <${ABSENCES_FILE} > ${ABSENCES_FILE}.fixed

if [ -d "${PROJECT_DIR}/pav_validation" ]; then
    rm -r ${PROJECT_DIR}/pav_validation
fi

mkdir ${PROJECT_DIR}/pav_validation

cd ${PROJECT_DIR}/pav_validation

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

# validate missing ogs using tblastn and blastp
while read -r LOST_GENOME LOST_OG; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}.fa -g ${LOST_GENOME_DIR}/${LOST_GENOME} \
        -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -d ${BLAST_DB} >> jobqueue
done < ${ABSENCES_FILE}.fixed

# parallelization stuff
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

## after everything is done, remove verbose outputs
awk -v OFS='\t' '{ if ($3 == "yes") {print $1, $2, $3} else {print $1, $2, "no"}}' ${OUTPUT_FILE} > ${OUTPUT_FILE}.simplified