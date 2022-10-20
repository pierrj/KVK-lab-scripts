#!/bin/bash
#SBATCH --job-name=gnu_parallel_rf_low_prio
#SBATCH --partition=savio
#SBATCH --qos=savio_lowprio
#SBATCH --account=co_minium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd ${SLURM_SUBMIT_DIR}

module load parallel

source activate /global/scratch/users/pierrj/conda_envs/random_forest/

parallel -j ${SLURM_NTASKS} < jobqueue_${node} > ${OUTPUT_FILE}_${node}