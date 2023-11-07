#!/bin/bash
#SBATCH --job-name=space_sniffer
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/projects/vector_kvklab/pierrj/

# du -h eccdna_backup > du_output_5_26_23

# rm -r eccdna_backup.tar.gz

# tar -czvf eccdna_backup.tar.gz eccdna_backup

rm -r /global/scratch/projects/vector_kvklab/pierrj/eccdna_backup

rm -r /global/scratch/users/pierrj/eccDNA