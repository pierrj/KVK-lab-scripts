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
mkdir /global/scratch/users/pierrj/repeatmasker/moryzae_pyrettest/
cd /global/scratch/users/pierrj/repeatmasker/moryzae_pyrettest/
RepeatMasker -pa 24 -s -nolow -no_is -norna -lib /global/scratch/users/pierrj/references/pyret_3group_test.fasta -gff -a -dir . /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta