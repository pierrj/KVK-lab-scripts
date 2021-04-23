#!/bin/bash
#SBATCH --job-name=ccs_demux
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
cd ${SLURM_SUBMIT_DIR}
cd ${sample}

source activate python=2.7

ccs -j ${SLURM_NTASKS} --reportFile=${sample}.ccs.report ${sample}_subreads.bam ${sample}_ccs.bam
lima --different --dump-removed --split-bam-named --ccs ${sample}_ccs.bam ${barcodes} ${sample}_demux.bam