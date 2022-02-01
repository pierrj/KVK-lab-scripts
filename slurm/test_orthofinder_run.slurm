#!/bin/bash
#SBATCH --job-name=test_orthofinder_run
#SBATCH --partition=savio2_knl
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/orthofinder_runs/gladieux_test_run

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

# orthofinder -oa -M msa -A mafft -T fasttree -S diamond_ultra_sens -t ${SLURM_NTASKS} -a 16 -o orthofinder_1_28_22 -f proteomes_combined_for_orthofinder_corrected

# orthofinder -op -S diamond_ultra_sens -f proteomes_combined_for_orthofinder_corrected -o orthofinder_1_28_22 > jobqueue

diamond blastp -d orthofinder_1_28_22/Results_Feb01/WorkingDirectory/diamond_ultra_sensDBSpecies55 -q orthofinder_1_28_22/Results_Feb01/WorkingDirectory/Species55.fa -o orthofinder_1_28_22/Results_Feb01/WorkingDirectory/Blast55_55.txt --ultra-sensitive -p 1 --quiet -e 0.001 --compress 1