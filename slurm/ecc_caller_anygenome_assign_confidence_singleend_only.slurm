#!/bin/bash
#SBATCH --job-name=doesnt_matter
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd ${SLURM_SUBMIT_DIR}
cd ${sample}

/global/home/users/pierrj/git/bash/assign_confidence_singleend_nomerge.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b no_secondary.filtered.sorted.${sample}.bam \
    -r ${sample}.confirmedsplitreads.bed