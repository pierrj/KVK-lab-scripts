#!/bin/bash
#SBATCH --job-name=repeatmasker
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/repeatmasker/ORSA_w_organelles
RepeatMasker -pa 24 -s -lib /global/scratch/users/pierrj/repbase/TIGR_Oryza_Repeats.v3.3 -gff -a -dir . /global/scratch/users/pierrj/references/ORSA_w_organelles.fasta