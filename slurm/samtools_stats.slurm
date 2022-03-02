#!/bin/bash
#SBATCH --job-name=samtools_stats
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_1

samtools stats SMRTcell_1_subreads.bam

samtools stats SMRTcell_1_ccs.bam

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_2

samtools stats SMRTcell_2_subreads.bam

samtools stats SMRTcell_2_ccs.bam

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_3

samtools stats SMRTcell_3_subreads.bam

samtools stats SMRTcell_3_ccs.bam