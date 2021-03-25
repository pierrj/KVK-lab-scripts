#!/bin/bash
#SBATCH --job-name=snakemake-sv-callers
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv-callers/snakemake

source ~/.bashrc
conda activate wf

snakemake -R -C echo_run=0 --use-conda --latency-wait 30 --jobs 14 --cluster "sbatch --partition=savio2 --account=fc_kvkallow --qos=savio_normal--nodes=1 --ntasks-per-node=24 --cpus-per-task={threads} --time=72:00:00 --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out --job-name=smk.{rule} --mail-user=pierrj@berkeley.edu --mail-type=ALL"