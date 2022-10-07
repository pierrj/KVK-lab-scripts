#!/bin/bash
#SBATCH --job-name=htc4_gnu_parallel_pannzer
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd ${SLURM_SUBMIT_DIR}

source activate /global/scratch/users/pierrj/conda_envs/pannzer

parallel -j ${SLURM_NTASKS} < jobqueue_${node}