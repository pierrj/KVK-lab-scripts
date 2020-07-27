#!/bin/bash
#SBATCH --job-name=normalize_and_average_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t .genomecoverage.filtered.bed -n filtered.sorted.

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalize_and_average_tests

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/mapfile_for_normalize_and_average -f 1000000 -b 100