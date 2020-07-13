#!/bin/bash
#SBATCH --job-name=lists_confirm_memory_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/confirmsrs_numpy_test
THREADS='20'
SAMPLE='CT_1A'
chrom_count='56'

# sort -k4,4 -k2,2n outwardfacing.CT_1A.renamed.bed > sorted.outwardfacing.CT_1A.renamed.bed


# awk -v OFS='\t' '{
#     prev=$0; f2=$2; f3=$3; f4=$4
#     getline 
#     if ($4 == f4 && f2 < $2 && f3 <$3) {
#         print $1, f2, $3, f3, $2
#     }
# }' sorted.outwardfacing.CT_1A.renamed.bed > sorted.grouped.outwardfacing.CT_1A.renamed.bed

ipcluster start -n ${THREADS} --cluster-id="cluster-id-${SAMPLE}" &
sleep 60
ipython /global/home/users/pierrj/git/python/ecc_caller_anygenome_confirmsrs.py lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} ${chrom_count}
ipcluster stop --cluster-id="cluster-id-${SAMPLE}"