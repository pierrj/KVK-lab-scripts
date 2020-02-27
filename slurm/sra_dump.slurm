#!/bin/bash
#SBATCH --job-name=sra_fastqdump
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/2018_moller
while read p; do
  /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump -O /global/scratch/users/pierrj/sra "$p"
done < SRR_list.txt
fastqc *.fastq --outdir=fastqc