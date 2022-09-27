#!/bin/bash
#SBATCH --job-name=file_cleanup
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap

rm -r predicted_effectors

rm -r *.effectorpout
rm -r *.notm.faa
rm -r *.notm.names
rm -r *.notm.sp.faa
rm -r *.predicted_effectors
rm -r *.signalpout.short
rm -r *.sp.names
rm -r *.tmhmmout
rm -r TMHMM_*