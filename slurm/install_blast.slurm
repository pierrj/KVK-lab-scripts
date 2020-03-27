#!/bin/bash
#SBATCH --job-name=download_blast
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/blast/
tar xvfz ncbi-blast-2.10.0+-src.tar.gaz
cd ncbi-blast-2.10.0+-src
cd c++
./configure
cd ReleaseMT/build
make all_r