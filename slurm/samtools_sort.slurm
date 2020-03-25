#!/bin/bash
#SBATCH --job-name=samtools_sort
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/guy11_and_spikeins
samtools sort mergedandpe.PQ_1A_bwamem.bam > sorted.mergedandpe.PQ_1A_bwamem.bam
samtools index sorted.mergedandpe.PQ_1A_bwamem.bam