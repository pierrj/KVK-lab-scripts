#!/bin/bash
#SBATCH --job-name=orthofinder_blast
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

cd ${SLURM_SUBMIT_DIR}

module purge

export PERL5LIB=''

parallel -j ${SLURM_NTASKS} < jobqueue_${node}