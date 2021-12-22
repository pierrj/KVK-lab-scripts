#!/bin/bash
#SBATCH --job-name=tblastn_validation_optimization
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/exonerate_validation/false_positive_rate

LOST_GENOME=/global/scratch/users/pierrj/moryzae_pav/exonerate_validation/CH0043.fasta
E_VALUE=1e-10
HIT_COUNT=2

if [ -f "optimization_out" ]; then
    rm optimization_out
fi

shuf og_list_all_presences_CH0043 | head -10000 > test_10000

for QUERY_COV in {60..70..5}
do
  for PIDENT in {60..80..5}
  do

    if [ -f "jobqueue" ]; then
        rm jobqueue
    fi

    while read og; do
        echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l orthogroup_protein_out_no_CH0043/${og}_protein.fasta -g ${LOST_GENOME} -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} >> jobqueue
    done < test_1000

    echo ${QUERY_COV}_${PIDENT}

    if [ -f "${QUERY_COV}_${PIDENT}_out" ]; then
        rm ${QUERY_COV}_${PIDENT}_out
    fi


    parallel -j 24 < jobqueue >> ${QUERY_COV}_${PIDENT}_out

    yes_count=$(awk '$2 == "yes"' ${QUERY_COV}_${PIDENT}_out | wc -l)
    no_count=$(awk '$2 == "no"' ${QUERY_COV}_${PIDENT}_out | wc -l)

    echo -e ${QUERY_COV}'\t'${PIDENT}'\t'$yes_count'\t'$no_count >> optimization_out

  done
done