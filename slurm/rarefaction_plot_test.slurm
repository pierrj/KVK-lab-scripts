#!/bin/bash
#SBATCH --job-name=rarefaction_plot_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/rarefaction_test

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
sample="RC_1A"
bamfile="/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/filtered.sorted.RC_1A.bam"

/global/home/users/pierrj/git/bash/generate_rarefaction_plot_data.sh -m ${mapfile} -s ${sample} -t ${SLURM_NTASKS} -b ${bamfile}