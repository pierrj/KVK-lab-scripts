#!/bin/bash
#SBATCH --job-name=anygenome_coverageplot
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/anygenome_coverageplot_test
/global/home/users/pierrj/git/bash/anygenome_coverageplot.sh -g /global/scratch/users/pierrj/references/Magnaporthe_oryzae.MG8.dna.toplevel.noMT_bwa -m map \
    -1 /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R1.fastq -2 /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R2.fastq \
    -s CT_1A -t 24