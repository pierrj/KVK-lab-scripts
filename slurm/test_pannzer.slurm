#!/bin/bash
#SBATCH --job-name=test_pannzer
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


/global/scratch/users/pierrj/conda_envs/pannzer/bin/python \
/global/scratch/users/pierrj/pannzer2/SANSPANZ.3/runsanspanz.py \
-R -o ",DE.out,GO.out,anno.out" -s "Pyricularia oryzae" < Magnaporthe_oryzae.MG8.pep.all.fa