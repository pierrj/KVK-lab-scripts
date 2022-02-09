#!/bin/bash
#SBATCH --job-name=fasttree_scos
#SBATCH --partition=savio_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/PAV_SV/PAV/full_run_tree_1_31_2022

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

SCO_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022/orthofinder_1_31_2022/Results_Feb03/WorkingDirectory/OrthoFinder/Results_Feb03/Single_Copy_Orthologue_Sequences
basename=$(basename ${SCO_DIR})
MSA=${basename}_msas_combined.trimmed.fasta

fasttree -gamma -out ${basename}_msas_combined.trimmed.fasttree ${basename}_msas_combined.trimmed.fasta