#!/bin/bash
#SBATCH --job-name=illumina_v_pacbio
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/illumina_v_pacbio

sample=G3_1A
mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"

/global/home/users/pierrj/git/bash/ecc_caller_mapq0_nomatchlength.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b uniq.filtered.sorted.${sample}.bam \
    -q multimapped.filtered.name_sorted.${sample}.bam

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b no_secondary.filtered.sorted.${sample}.bam \
    -r ${sample}.confirmedsplitreads.bed