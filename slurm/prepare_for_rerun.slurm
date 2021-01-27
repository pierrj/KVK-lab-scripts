#!/bin/bash
#SBATCH --job-name=prepare_for_rerun
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

while read sample; do
cd ${sample}
    mv ${sample}_R1.fastq ../${sample}_R1.fastq
    mv ${sample}_R2.fastq ../${sample}_R1.fastq
cd ..
done < mapfile

while read sample; do
cd ${sample}
    rm -r *
cd ..
done < mapfile


while read sample; do
cd ${sample}
    mv ../${sample}_R1.fastq ${sample}_R1.fastq
    mv ../${sample}_R2.fastq ${sample}_R2.fastq
cd ..
done < mapfile