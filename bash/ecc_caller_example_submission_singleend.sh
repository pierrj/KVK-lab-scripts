#!/bin/bash
mapfile="/global/scratch/users/pierrj/references/Scer_S288C.contignames"
genome_bwa="/global/scratch/users/pierrj/references/Scer_S288C_bwa"
sample="ERR525229"
READONE_LINK="ftp.sra.ebi.ac.uk/vol1/fastq/ERR525/ERR525229/ERR525229.fastq.gz"
sbatch --job-name=$sample.ecc_caller_single_end --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK /global/home/users/pierrj/git/slurm/ecc_caller_singleend.slurm