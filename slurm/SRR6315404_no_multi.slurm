#!/bin/bash
#SBATCH --job-name=SRR6315404_no_multi
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human

sample=SRR6315404
mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"

/global/home/users/pierrj/git/bash/ecc_caller_only_unique.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/uniq.filtered.sorted.${sample}.bam

/global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/no_secondary.filtered.sorted.${sample}.bam \
    -r ${sample}.confirmedsplitreads.bed