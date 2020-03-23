#!/bin/bash
#SBATCH --job-name=bwaindex
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
cd /global/scratch/users/pierrj/references
bwa index -p guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_bwa guy11_genome_baoetal2017_with_70-15_mito_with_spikeins.fasta