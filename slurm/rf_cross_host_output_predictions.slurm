#!/bin/bash
#SBATCH --job-name=rf_cross_host_output_predictions
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/cross_host

source activate /global/scratch/users/pierrj/conda_envs/random_forest/

MAJORITY_FRACTION=0.5
APPROACH=RF
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=None
DEPTH=None
BOOTSTRAP=True

INPUT_DF=gene_info.cross_host.wheat_blast.txt
INPUT_DF2=gene_info.cross_host.rice_blast.txt
OUTPUT_STRING=cross_host_wheat_first
OUTPUT_FILE=rf_results.${OUTPUT_STRING}.txt

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_cross_host_output_preds.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP $INPUT_DF2 $OUTPUT_STRING > $OUTPUT_FILE


INPUT_DF=gene_info.cross_host.rice_blast.txt
INPUT_DF2=gene_info.cross_host.wheat_blast.txt
OUTPUT_STRING=cross_host_rice_first
OUTPUT_FILE=rf_results.${OUTPUT_STRING}.txt

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_cross_host_output_preds.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP $INPUT_DF2 $OUTPUT_STRING > $OUTPUT_FILE