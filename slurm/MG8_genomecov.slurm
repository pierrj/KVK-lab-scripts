#!/bin/bash
#SBATCH --job-name=ecc_caller_MG8_genomecov
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/CT_1A_MG8
sample=CT_1A
bedtools genomecov -d -i parallel.confirmed.sorted.bed -g /global/scratch/users/pierrj/references/Magnaporthe_oryzae.MG8.dna.toplevel.noMT.chromsizes > parallel.confirmed.sorted.bed.genomecov
awk -v OFS='\t' '{sum+=$3} NR%100==0 {print $1, $2, sum/100; sum =0}'  parallel.confirmed.sorted.bed.genomecov > parallel.confirmed.sorted.bed.manhattanplot.old
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT' > first_row
cat first_row parallel.confirmed.sorted.bed.manhattanplot.old > parallel.confirmed.sorted.bed.manhattanplot