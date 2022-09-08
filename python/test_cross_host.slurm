#!/bin/bash
#SBATCH --job-name=test_cross_host
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

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_cross_host_test.py \
                                    /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/random_forest/gene_info_w_proportions.csv \
                                    /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/random_forest/gene_info_w_proportions.csv \
                                    train_rice_test_wheat > train_rice_test_wheat.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/random_forest

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_cross_host_test.py \
                                    /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/random_forest/gene_info_w_proportions.csv \
                                    /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/random_forest/gene_info_w_proportions.csv \
                                    train_wheat_test_rice > train_wheat_test_rice.out