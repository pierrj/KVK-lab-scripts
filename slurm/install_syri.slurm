#!/bin/bash
#SBATCH --job-name=install_syri
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

export CONDA_PKGS_DIRS=/global/scratch/users/pierrj/conda_envs/

source activate /global/scratch/users/pierrj/conda_envs/syri_2

conda install -y cython numpy scipy pandas=0.23.4 biopython psutil matplotlib=3.0.0
conda install -y -c conda-forge python-igraph
conda install -y -c bioconda pysam