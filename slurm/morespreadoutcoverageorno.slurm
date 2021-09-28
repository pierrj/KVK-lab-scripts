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


# cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A

# GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

# bedtools makewindows -w 100 -g ${GENOME_FILE} > guy11.genome_100bp_windows.bed

# bedtools coverage -sorted -a guy11.genome_100bp_windows.bed -b no_secondary.filtered.sorted.G3_1A.bam > guy11.genome_100bp_windows.coverage


# awk '{sum += $4} END { print sum/NR }' guy11.genome_100bp_windows.coverage

# 468.152

# awk '{x+=$4;y+=$4^2}END{print sqrt(y/NR-(x/NR)^2)}' guy11.genome_100bp_windows.coverage

# 3922.62



GENOME_FILE=/global/scratch/users/pierrj/references/GRCh37.p13.chromsizes
SAMPLE=SRR6315399
ORG=human

cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${SAMPLE}

bedtools makewindows -w 100 -g ${GENOME_FILE} > ${ORG}.genome_100bp_windows.bed

bedtools coverage -sorted -a ${ORG}.genome_100bp_windows.bed -b no_secondary.filtered.sorted.${SAMPLE}.bam > ${ORG}.genome_100bp_windows.coverage


awk '{sum += $4} END { print sum/NR }' ${ORG}.genome_100bp_windows.coverage

awk '{x+=$4;y+=$4^2}END{print sqrt(y/NR-(x/NR)^2)}' ${ORG}.genome_100bp_windows.coverage