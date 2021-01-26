#!/bin/bash
#SBATCH --job-name=ecc_calling_mapq0_implementation_fulltest
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_with_mapq0/implementation_fulltest

GENOME_DB=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa
READONE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R1.fastq
READTWO=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R2.fastq
SAMPLE=G3_1A
MAPFILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
FILTERED_BAMFILE=uniq.filtered.sorted.G3_1A.bam
FILTERED_BAMFILE_QSORTED=multimapped.filtered.name_sorted.G3_1A.bam
CONFIRMED_SPLITREAD=G3_1A.confirmedsplitreads.bed
NO_SECONDARY=no_secondary.filtered.sorted.G3_1A.bam

/global/home/users/pierrj/git/bash/generate_bam_file_mapq0.sh -g ${GENOME_DB} \
    -1 ${READONE} \
    -2 ${READTWO} \
    -s ${SAMPLE} \
    -t ${SLURM_NTASKS} \
    -m ${MAPFILE}

/global/home/users/pierrj/git/bash/ecc_caller_mapq0.sh -m ${MAPFILE} \
    -s ${SAMPLE} \
    -t ${SLURM_NTASKS} \
    -b ${FILTERED_BAMFILE} \
    -q ${FILTERED_BAMFILE_QSORTED}

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${MAPFILE} \
    -s ${SAMPLE} \
    -t ${SLURM_NTASKS} \
    -b ${NO_SECONDARY} \
    -r ${CONFIRMED_SPLITREAD}