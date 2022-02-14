#!/bin/bash
#SBATCH --job-name=raxml_bootstrapping
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/raxml_ng_test

/global/scratch/users/pierrj/raxml_ng_savio1/bin/raxml-ng-mpi --bootstrap --msa savio1_T1.raxml.rba --prefix BS1 --seed 1111 --threads 10 --bs-trees 5