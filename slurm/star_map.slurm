#!/bin/bash
#SBATCH --job-name=star_map
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
cd /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes/47_star
STAR --runThreadN 24 \
    --genomeDir /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed_star \
    --readFilesIn /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes/SRR298847.sra.fastq \
    --quantMode GeneCounts
