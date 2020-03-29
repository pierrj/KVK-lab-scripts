#!/bin/bash
#SBATCH --job-name=untar_everything
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
cd /global/scratch/users/pierrj/blast/blastdb_ver4_curl/blast/db/v4
# wget -r -l1 --no-parent --no-check-certificate -A nt_v4.*.tar.gz* 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/v4/' -P /global/scratch/users/pierrj/blast/blastdb_ver4_curl
for file in *.tar.gz; do
    tar xvfz $file
    size=${#file}
    mv ${file:0:size-7}/* .
    rm ${file:0:size-7}
done