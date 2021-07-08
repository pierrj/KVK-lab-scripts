#!/bin/bash
#SBATCH --job-name=length_dsn
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/length_dsn

sample=G3_1A
mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"


/global/home/users/pierrj/git/bash/ecc_caller_mapq0_no_rm.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/uniq.filtered.sorted.${sample}.bam \
    -q /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/multimapped.filtered.name_sorted.${sample}.bam

awk '{print $3-$2}' unique_parallel.confirmed > unique_parallel.confirmed.length_dsn
awk '{print $3-$2}' mapq0_parallel.confirmed > mapq0_parallel.confirmed.length_dsn
