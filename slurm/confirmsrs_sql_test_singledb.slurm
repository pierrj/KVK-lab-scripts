#!/bin/bash
#SBATCH --job-name=confirmsrs_sql_test_singledb
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sqlconfirmsrs_rice
THREADS='20'
SAMPLE='RC_1A'
chrom_count='57'

ipcluster start -n ${THREADS} --cluster-id="cluster-id-${SAMPLE}" &
sleep 300
ipython /global/home/users/pierrj/git/python/ecc_caller_anygenome_confirmsrs_sql_singledb.py lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} ${chrom_count}
ipcluster stop --cluster-id="cluster-id-${SAMPLE}"