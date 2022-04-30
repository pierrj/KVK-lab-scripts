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


cd /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/optimization_run

LOST_GENOME=/global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/CH0043.fasta
E_VALUE=1e-10
HIT_COUNT=2
QUERY_COV=55
PIDENT=55
GENE_GFF=/global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate_fungap_or_orthofinder/CH0043_merge.gff3

if [ -f "optimization_out" ]; then
    rm optimization_out
fi

# for GENE_OVERLAP in {80..95..5}
# do

GENE_OVERLAP=100
if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read og; do
    echo /global/home/users/pierrj/git/bash/tblastn_validation.sh -l /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate/orthogroup_protein_out_no_CH0043/${og}_protein.fasta \
    -g ${LOST_GENOME} -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -f  ${GENE_GFF} -b ${GENE_OVERLAP} >> jobqueue
done < /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/false_positive_rate/og_list_all_presences_CH0043

echo ${GENE_OVERLAP}

if [ -f "${GENE_OVERLAP}_out" ]; then
    rm ${GENE_OVERLAP}_out
fi


parallel -j 24 < jobqueue >> ${GENE_OVERLAP}_out

yes_otherog_count=$(awk '$3 == "yes_otherog"' ${GENE_OVERLAP}_out | wc -l)
yes_unannotated_count=$(awk '$3 == "yes_unannotated"' ${GENE_OVERLAP}_out | wc -l)
no_count=$(awk '$3 == "no"' ${GENE_OVERLAP}_out | wc -l)

echo -e ${GENE_OVERLAP}'\t'$yes_otherog_count'\t'$yes_unannotated_count'\t'$no_count >> optimization_out
# done