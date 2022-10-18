#!/bin/bash
#SBATCH --job-name=htc4_gnu_parallel_rf_single_cpu
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --account=co_minium
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd ${SLURM_SUBMIT_DIR}

module load parallel

source activate /global/scratch/users/pierrj/conda_envs/random_forest/

parallel -j ${SLURM_NTASKS} < jobqueue_${cpu} >> ${OUTPUT_FILE}