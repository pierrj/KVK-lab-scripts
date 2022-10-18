#!/bin/bash
#SBATCH --job-name=rf_cross_host
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks-per-node=2
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/cross_host

MAJORITY_FRACTION=1.0
APPROACH=SMOTE
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=sqrt
DEPTH=None
BOOTSTRAP=True

INPUT_DF=gene_info.rice_blast.txt
INPUT_DF2=gene_info.wheat_blast.txt
OUTPUT_FILE=rf_results.cross_host_rice_first.txt

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP $INPUT_DF2


INPUT_DF=gene_info.wheat_blast.txt
INPUT_DF2=gene_info.rice_blast.txt
OUTPUT_FILE=rf_results.cross_host_wheat_first.txt

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP $INPUT_DF2