#!/bin/bash
#SBATCH --job-name=filter_mapq0
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
while read sample; do
samtools view -b -q 1 ${sample}/${sample}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_with_mapq0/magnaporthe_pure_culture/${sample}/filtered.sorted.${sample}.bam
done < mapfile