#!/bin/bash
#SBATCH --job-name=ecc_density_counting
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

OUTPUT=G3_ecc_density.txt
GENOME=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
MAPFILE=mapfile

/global/home/users/pierrj/git/bash/ecc_density_counting.sh -o ${OUTPUT} -g ${GENOME} -m ${MAPFILE}