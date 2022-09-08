#!/bin/bash
#SBATCH --job-name=test_perf_drop_genomes
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/random_forest

module purge

source activate /global/scratch/users/pierrj/conda_envs/random_forest/

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_test_perf_drop_genomes.py gene_info_w_proportions.csv test_perf_drop_genomes > test_perf_drop_genomes.out