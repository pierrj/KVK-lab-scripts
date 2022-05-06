#!/bin/bash
#SBATCH --job-name=gnu_parallel_multinode
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd ${SLURM_SUBMIT_DIR}

source activate /global/scratch/users/pierrj/conda_envs/agat
PERL5LIB=''

parallel -j ${SLURM_NTASKS} < jobqueue_${node} >> ${OUTPUT_FILE}