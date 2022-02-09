#!/bin/bash
#SBATCH --job-name=raxml_multinode
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/full_run_tree_1_31_2022

module purge
module load raxml
module load gcc openmpi

SCO_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022/orthofinder_1_31_2022/Results_Feb03/WorkingDirectory/OrthoFinder/Results_Feb03/Single_Copy_Orthologue_Sequences
basename=$(basename ${SCO_DIR})
MSA=${basename}_msas_combined.trimmed.fasta
output=gladieux_full_raxml_multinode

mpirun raxmlHPC-HYBRID-SSE3 -s ${MSA} -n ${output} -m PROTGAMMAGTR -T 10 -f a -x 12345 -p 12345 -# 100