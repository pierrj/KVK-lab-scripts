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

ls orthogroup_protein_out_no_CH0043 | awk '{print substr($1, 0, 9)}' > no_CH0043_there

shuf og_list_all_presences_CH0043 | head -1000 > test_1000

cat no_CH0043_there test_1000 | sort | uniq -d > test_a_few


# if [ -f "jobqueue" ]; then
#     rm jobqueue
# fi

# while read og; do
#     echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l orthogroup_protein_out_no_CH0043/${og}_protein.fasta -g ${LOST_GENOME} -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} >> jobqueue
# done < test_a_few

# parallel -j 24 < jobqueue


while read og; do
    /global/home/users/pierrj/git/bash/tblastn_validation.sh -l orthogroup_protein_out_no_CH0043/${og}_protein.fasta -g ${LOST_GENOME} -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT}
done < test_a_few