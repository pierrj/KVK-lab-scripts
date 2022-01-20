#!/bin/bash
#SBATCH --job-name=install_fungap
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

module purge

source activate /global/scratch/users/pierrj/conda_envs/fungap

conda install -y braker2=2.1.5 trinity=2.12.0 repeatmodeler=2.0.1 hisat2=2.2.1 pfam_scan=1.6 busco=5.1.2

pip install biopython bcbio-gff markdown2 matplotlib

cpanm YAML Hash::Merge Logger::Simple Parallel::ForkManager MCE::Mutex Thread::Queue threads