#!/bin/bash
#SBATCH --job-name=mapq0_ecc_caller_2018_moller
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/mapq0_2018_moller

mapfile='/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames'
sample='SRR6315407'

/global/home/users/pierrj/git/bash/call_ecc_regions_allmapq.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b /global/scratch/users/pierrj/eccDNA/2018_moller/${sample}/${sample}.sorted.mergedandpe.bwamem.bam