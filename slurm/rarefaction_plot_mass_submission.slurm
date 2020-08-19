#!/bin/bash
#SBATCH --job-name=doesnt_matter
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
cd ${SLURM_SUBMIT_DIR}
cd ${sample}

bamfile=${sample}.sorted.mergedandpe.bwamem.bam

/global/home/users/pierrj/git/bash/generate_rarefaction_plot_data.sh -m ${mapfile} -s ${sample} -t ${SLURM_NTASKS} -b ${bamfile}