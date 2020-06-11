#!/bin/bash
#SBATCH --job-name=copia_srs
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/copia_srs/
sample=G3_1A

bedtools intersect -f 1 -wa -a SRs.mergedandpe.${sample}_bwamem.bed -b /global/scratch/users/pierrj/references/GUY11_pacbio_copia_group1_ltr_locs.bed > copia.SRs.mergedandpe.${sample}_bwamem.bed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' copia.SRs.mergedandpe.${sample}_bwamem.bed copia.SRs.mergedandpe.${sample}_bwamem.bed > exactlytwice.copia.SRs.mergedandpe.${sample}_bwamem.bed