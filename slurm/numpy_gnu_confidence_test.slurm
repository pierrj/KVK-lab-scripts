#!/bin/bash
#SBATCH --job-name=confidence_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/numpy_gnu_confidence_test

/global/home/users/pierrj/git/bash/assign_confidence.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames\
    -s CT_1A -t 20 -c CT_1A.genomecoverage.bed

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sql_confidence/slow

THREADS='20'
SAMPLE='CT_1A_test'
chrom_count='56'

ipcluster start -n ${THREADS} --cluster-id="cluster-id-${SAMPLE}" &
sleep 60
ipython /global/home/users/pierrj/git/python/ecc_caller_anygenome_assignconfidence_slow.py ${SAMPLE}.genomecoverage.filtered.renamed.bed ${SAMPLE} 56
ipcluster stop --cluster-id="cluster-id-${SAMPLE}"