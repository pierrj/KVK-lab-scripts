#!/bin/bash
#SBATCH --job-name=junction_split_reads_per_sample
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
        echo ${sample}
        samtools view -c no_secondary.filtered.sorted.${sample}.bam
        wc -l ${sample}.ecc_caller_out.splitreads.bed
    cd ..
done < mapfile