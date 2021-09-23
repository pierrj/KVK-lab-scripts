#!/bin/bash
#SBATCH --job-name=morespreadoutcoverageorno
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A

GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

bedtools makewindows -w 100 -g ${GENOME_FILE} > guy11.genome_100bp_windows.bed

bedtools coverage -a guy11.genome_100bp_windows.bed -b no_secondary.filtered.sorted.G3_1A.bam > guy11.genome_100bp_windows.coverage