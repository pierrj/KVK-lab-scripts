#!/bin/bash
#SBATCH --job-name=ltr_srs_covfile
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/LTR_splitreads
bedtools genomecov -d -ibam SRs.no_orientation.allmapq.OsSR.bam > genomecov.SRs.no_orientation.allmapq.OsSR.bam