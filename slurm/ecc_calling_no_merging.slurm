#!/bin/bash
#SBATCH --job-name=ecc_calling_no_merging
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_no_merging

THREADS=${SLURM_NTASKS}
SAMPLE=G3_1A

parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_nodb_nomerge.py ${SAMPLE} {} renamed.filtered.sorted.${SAMPLE}.bam ::: $(seq -w 1 ${THREADS})