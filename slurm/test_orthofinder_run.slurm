#!/bin/bash
#SBATCH --job-name=test_orthofinder_run
#SBATCH --partition=savio2_knl
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/orthofinder_runs/gladieux_test_run

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

orthofinder -oa -M msa -A mafft -T fasttree -S diamond_ultra_sens -t ${SLURM_NTASKS} -a 16 -o orthofinder_1_28_22 -f proteomes_combined_for_orthofinder_corrected