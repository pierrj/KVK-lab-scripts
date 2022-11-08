#!/bin/bash
#SBATCH --job-name=dummy_job
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --cpus=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


echo hello