#!/bin/bash
#SBATCH --job-name=sco_tree_making_optimization
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/tree_making_optimization

module load raxml

date +"%T"

raxmlHPC-HYBRID-SSE3 -s ALL_SCOs.trim.afa -n raxmlHPC-HYBRID-SSE3.GAMMA.ALL_SCOs -m PROTGAMMAGTR -T 24 -f a -x 12345 -p 12345 -# 100

date +"%T"

date +"%T"

raxmlHPC-HYBRID-SSE3 -s ALL_SCOs.trim.afa -n raxmlHPC-HYBRID-SSE3.CAT.ALL_SCOs -m PROTCATGTR -T 24 -f a -x 12345 -p 12345 -# 100

date +"%T"