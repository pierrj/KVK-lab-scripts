#!/bin/bash
#SBATCH --job-name=bedtools_genome_cov
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11

while read sample; do
    bedtools genomecov -d -ibam ${SRA}.Aligned.sortedByCoord.out.bam > ${SRA}.genomecov
done < zhang_et_al_2019_complete_medium

while read sample; do
    bedtools genomecov -d -ibam ${SRA}.Aligned.sortedByCoord.out.bam > ${SRA}.genomecov
done < zhang_et_al_2019_in_planta

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/b71

while read sample; do
    bedtools genomecov -d -ibam ${SRA}.Aligned.sortedByCoord.out.bam > ${SRA}.genomecov
done < b71_assembly_complete_medium

while read sample; do
    bedtools genomecov -d -ibam ${SRA}.Aligned.sortedByCoord.out.bam > ${SRA}.genomecov
done < b71_assembly_complete_medium_in_planta