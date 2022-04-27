#!/bin/bash
#SBATCH --job-name=process_fungap_out
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes


if [ ! -d "all_proteomes_corrected" ]; then
    mkdir all_proteomes_corrected
fi

source activate /global/scratch/users/pierrj/conda_envs/orthofinder
module purge

python /global/home/users/pierrj/git/python/process_protein_sequences_for_orthofinder.py all_proteomes \
    all_proteomes_corrected \
    gladieux_et_al_assembly_stats.txt \
    output_accession_names