#!/bin/bash
#SBATCH --job-name=wgs_count_comparison
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/wgs_data/


mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample="ERR2660591"
sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample="ERR5060418"
sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample="SRR11528297"
sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm

awk '$5 != "lowq"' ERR2660591/ERR2660591.ecc_caller_out.details.txt | wc -l

wc -l ERR2660591/ERR2660591_sra_1.fastq | awk '{print $1*2}'

awk '$5 != "lowq"' ERR5060418/ERR5060418.ecc_caller_out.details.txt | wc -l

wc -l ERR5060418/ERR5060418_sra_1.fastq | awk '{print $1*2}'

awk '$5 != "lowq"' SRR11528297/SRR11528297.ecc_caller_out.details.txt | wc -l

wc -l SRR11528297/SRR11528297_sra_1.fastq | awk '{print $1*2}'