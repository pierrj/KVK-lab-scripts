#!/bin/bash
mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm; done < mapfile


mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample="ERR1830502"
READONE_LINK="ftp.sra.ebi.ac.uk/vol1/fastq/ERR183/002/ERR1830502/ERR1830502_1.fastq.gz"
READTWO_LINK="ftp.sra.ebi.ac.uk/vol1/fastq/ERR183/002/ERR1830502/ERR1830502_2.fastq.gz"
sbatch --job-name=$sample.ecc_caller_ena --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK,READTWO_LINK=$READTWO_LINK /global/home/users/pierrj/git/slurm/ecc_caller_ena.slurm

while IFS=$"\t" read -r sample READONE_LINK READTWO_LINK
do
echo ${sample}
echo $READONE_LINK
echo $READTWO_LINK
done < mapfile_ENA_ory

while IFS=$'\t' read -r sample READONE_LINK READTWO_LINK; do echo ${sample}; echo $READONE_LINK; echo $READTWO_LINK; done < mapfile_ENA_ory

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
while IFS=$'\t' read -r sample READONE_LINK READTWO_LINK; do sbatch --job-name=$sample.ecc_caller_ena --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK,READTWO_LINK=$READTWO_LINK /global/home/users/pierrj/git/slurm/ecc_caller_ena.slurm; done < mapfile_ENA_ory

mapfile="/global/scratch/users/pierrj/references/TAIR10.contignames"
genome_bwa="/global/scratch/users/pierrj/references/TAIR10_bwa"
while IFS=$'\t' read -r sample READONE_LINK READTWO_LINK; do sbatch --job-name=$sample.ecc_caller_ena --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK,READTWO_LINK=$READTWO_LINK /global/home/users/pierrj/git/slurm/ecc_caller_ena.slurm; done < mapfile_ENA_ara

mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
while IFS=$'\t' read -r sample READONE_LINK READTWO_LINK; do sbatch --job-name=$sample.ecc_caller_ena --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK,READTWO_LINK=$READTWO_LINK /global/home/users/pierrj/git/slurm/ecc_caller_ena.slurm; done < mapfile_ENA_ory


## post no merge


mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_assign_confidence_only.slurm; done < mapfile_ory

mapfile="/global/scratch/users/pierrj/references/TAIR10.contignames"
genome_bwa="/global/scratch/users/pierrj/references/TAIR10_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_assign_confidence_only.slurm; done < mapfile_ara