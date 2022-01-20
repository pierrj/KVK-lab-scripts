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

conda install -y python=3.8 cython=0.29.23 numpy=1.20.2 scipy=1.6.2 pandas=1.2.4 psutil=5.8.0 matplotlib=3.3.4
conda install -y python-igraph=0.9.1
conda install -y pysam=0.16.0.1
conda install -y plotsr

module purge

cd ~
git clone https://github.com/schneebergerlab/syri.git
cd syri
