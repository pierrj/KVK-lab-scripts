#!/bin/bash
#SBATCH --job-name=raxml_ng_scos
#SBATCH --partition=savio2_knl
#SBATCH --qos=savio_normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

## to compile 
#  wget https://github.com/amkozlov/raxml-ng/releases/download/1.1.0/raxml-ng_v1.1.0_linux_x86_64_MPI.zip
# unzip raxml-ng_v1.1.0_linux_x86_64_MPI.zip
# ./install.sh

module purge
module load gcc/7.4.0
module load openmpi
module load cmake

cd /global/scratch/users/pierrj/PAV_SV/PAV/raxml_ng_test

# /global/scratch/users/pierrj/raxml_ng/bin/raxml-ng-mpi --parse --msa Single_Copy_Orthologue_Sequences_msas_combined.trimmed.fasta --model PROTGTR+G --prefix T4

/global/scratch/users/pierrj/raxml_ng/bin/raxml-ng-mpi --msa T4.raxml.rba --prefix H2 --threads 64 --model PROTGTR+G