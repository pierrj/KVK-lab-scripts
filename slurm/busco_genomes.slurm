#!/bin/bash
#SBATCH --job-name=busco_genomes
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/references/moryzae_genomes_10_26_21/

source activate /global/scratch/users/pierrj/conda_envs/BUSCO

module remove python

busco -i fastas -l sordariomycetes_odb10 -o fastas_buscofied -m genome -c 24 -f