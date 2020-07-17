#!/bin/bash
mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11bao_7015mito_spikeins_tweaked_osativa_worganelles_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile



mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample="SRR6315407"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample="OsCT1"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm