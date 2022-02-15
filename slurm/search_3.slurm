#!/bin/bash
#SBATCH --job-name=search_3
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

module purge
module load gcc/7.4.0
module load openmpi
module load cmake

cd /global/scratch/users/pierrj/PAV_SV/PAV/raxml_ng_test

# /global/scratch/users/pierrj/raxml_ng_savio1/bin/raxml-ng-mpi --parse --msa Single_Copy_Orthologue_Sequences_msas_combined.trimmed.fasta --model PROTGTR+G --prefix savio1_T1

mpirun /global/scratch/users/pierrj/raxml_ng_savio1/bin/raxml-ng-mpi --msa savio1_T1.raxml.rba --prefix search_3 --threads 10 --extra thread-pin --seed 33333