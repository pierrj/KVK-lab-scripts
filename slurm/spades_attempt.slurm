#!/bin/bash
#SBATCH --job-name=spades_attempt
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/
module remove python
module load spades
spades.py -o spades --plasmid -1 CT_1A_R1.fastq -2 CT_1A_R2.fastq -m 512 -t 20