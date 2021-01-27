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
    rm -r *
cd ..
done < mapfile

while read sample; do
cd ${sample}
    cp ../../illumina_pre_mapq0_ecc_calling/${sample}/${sample}_R1.fastq .
    cp ../../illumina_pre_mapq0_ecc_calling/${sample}/${sample}_R2.fastq .
cd ..
done < mapfile