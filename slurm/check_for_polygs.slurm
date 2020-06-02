#!/bin/bash
#SBATCH --job-name=check_for_polyGs
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina
while read sample; 
do
    cd ${sample}
    echo ${sample}
    grep 'GGGGGGGGGGGGGGGGGGGG' ${sample}_R1.fastq | wc -l
    grep 'GGGGGGGGGGGGGGGGGGGG' ${sample}_R2.fastq | wc -l
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/
while read sample; 
do
    cd ${sample}
    echo ${sample}
    grep 'GGGGGGGGGGGGGGGGGGGG' ${sample}_R1.fastq | wc -l
    grep 'GGGGGGGGGGGGGGGGGGGG' ${sample}_R2.fastq | wc -l
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/
while read sample; 
do
    echo ${sample}
    grep 'GGGGGGGGGGGGGGGGGGGG' ${sample}_R1.fastq | wc -l
    grep 'GGGGGGGGGGGGGGGGGGGG' ${sample}_R2.fastq | wc -l
done < mapfile