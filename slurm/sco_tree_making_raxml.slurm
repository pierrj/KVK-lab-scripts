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

cd /global/scratch/users/pierrj/moryzae_pav/tree_making_optimization/mine

module load raxml
module load gcc openmpi

rm RAxML*

mpirun raxmlHPC-HYBRID-SSE3 -s SCO_nov_01_2021_combined_msas.trimmed.fasta -n SCO_nov_01_2021_combined_msas.trimmed.raxmltree -m PROTGAMMAGTR -T 10 -f a -x 12345 -p 12345 -# 100