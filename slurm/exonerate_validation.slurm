#!/bin/bash
#SBATCH --job-name=exonerate_validation
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/exonerate_validation

module load exonerate

LOST_GENOME=/global/scratch/users/pierrj/moryzae_pav/gladieux_et_al_2020_data/gladieux_et_al_2021_assemblies/CH0043.fasta
LOST_OG=/global/scratch/users/pierrj/moryzae_pav/orthogroup_cons/orthogroup_protein_out/OG0000323_protein.fasta

exonerate --model protein2genome ${LOST_OG} ${LOST_GENOME} --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > missing.gff