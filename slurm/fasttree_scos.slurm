#!/bin/bash
#SBATCH --job-name=fasttree_scos
#SBATCH --partition=savio_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

SCO_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/orthofinder_4_27_2022/Results_Apr28/WorkingDirectory/OrthoFinder/Results_Apr28/Single_Copy_Orthologue_Sequences
basename=$(basename ${SCO_DIR})
MSA=${basename}_msas_combined.trimmed.fasta

fasttree -gamma -out ${basename}_msas_combined.trimmed.fasttree ${basename}_msas_combined.trimmed.fasta