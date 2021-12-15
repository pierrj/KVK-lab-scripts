#!/bin/bash
#SBATCH --job-name=tblastn_validation_false_positive_tests
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
PIDENT=75
QUERY_COV=80
HIT_COUNT=2

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read og; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh orthogroup_protein_out_no_CH0043/${og} ${LOST_GENOME} ${E_VALUE} ${PIDENT} ${QUERY_COV} ${HIT_COUNT} >> jobqueue
done < og_list_all_presences_CH0043

parallel -j 24 < jobqueue