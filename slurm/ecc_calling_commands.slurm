#!/bin/bash
#SBATCH --job-name=ecc_calling_commands
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

# eccDNA calling from moryzae
mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm; done < mapfile


cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/

# eccDNA calling from osativa from this study
mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm; done < mapfile

cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run

# eccDNA calling from hsapiens data from moller et al 2018
mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm; done < mapfile

cd /global/scratch/users/pierrj/eccDNA/2020_wang/

# eccDNA calling from athaliana data from wang et al 2021
mapfile="/global/scratch/users/pierrj/references/TAIR10.contignames"
genome_bwa="/global/scratch/users/pierrj/references/TAIR10_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm; done < mapfile


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run/

# eccDNA calling from athaliana data from lanciano et al 2021
mapfile="/global/scratch/users/pierrj/references/TAIR10.contignames"
genome_bwa="/global/scratch/users/pierrj/references/TAIR10_bwa"
while IFS=$'\t' read -r sample READONE_LINK READTWO_LINK; do sbatch --job-name=$sample.ecc_caller_ena --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK,READTWO_LINK=$READTWO_LINK /global/home/users/pierrj/git/slurm/ecc_caller_ena.slurm; done < mapfile_ENA_ara

# eccDNA calling from osativa data from lanciano et al 2021
mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
while IFS=$'\t' read -r sample READONE_LINK READTWO_LINK; do sbatch --job-name=$sample.ecc_caller_ena --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa,READONE_LINK=$READONE_LINK,READTWO_LINK=$READTWO_LINK /global/home/users/pierrj/git/slurm/ecc_caller_ena.slurm; done < mapfile_ENA_ory

cd /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/

# eccDNA calling from scerevisiae data from moller et al 2015
mapfile="/global/scratch/users/pierrj/references/Scer_S288C.contignames"
genome_bwa="/global/scratch/users/pierrj/references/Scer_S288C_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller_single_end --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_singleend.slurm; done < mapfile