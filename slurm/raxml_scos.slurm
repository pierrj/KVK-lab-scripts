#!/bin/bash
#SBATCH --job-name=raxml_scos
#SBATCH --partition=savio2_knl
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/full_run_tree_1_31_2022

module load raxml

SCO_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022/orthofinder_1_31_2022/Results_Feb03/WorkingDirectory/OrthoFinder/Results_Feb03/Single_Copy_Orthologue_Sequences
basename=$(basename ${SCO_DIR})
MSA=${basename}_msas_combined.trimmed.fasta
output=gladieux_full_raxml

raxmlHPC-PTHREADS-SSE3 -s ${MSA} -n ${output} -m PROTGAMMAGTR -T ${SLURM_NTASKS} -f a -x 12345 -p 12345 -# 100