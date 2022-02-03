#!/bin/bash
#SBATCH --job-name=full_orthofinder_run
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022

orthofinder -oa -M msa -A mafft -T fasttree -t ${SLURM_NTASKS} -a ${SLURM_NTASKS} -b orthofinder_1_31_2022