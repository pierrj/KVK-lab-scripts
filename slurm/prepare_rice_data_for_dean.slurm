#!/bin/bash
#SBATCH --job-name=prepare_rice_data_for_dean
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
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/

while read sample; do
cd ${sample}
gzip -c ${sample}_R1.fastq > /global/home/groups/fc_kvkallow/${sample}_R1.fastq.gz
gzip -c ${sample}_R2.fastq > /global/home/groups/fc_kvkallow/${sample}_R2.fastq.gz
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/stress_experiments/infected/

while read sample; do
cd ${sample}
gzip -c ${sample}_R1.fastq > /global/home/groups/fc_kvkallow/${sample}_R1.fastq.gz
gzip -c ${sample}_R2.fastq > /global/home/groups/fc_kvkallow/${sample}_R2.fastq.gz
cd ..
done < mapfile