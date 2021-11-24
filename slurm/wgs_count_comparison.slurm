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



cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/wgs_data/

# call eccdnas with WGS data from previously publications and SRA
mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample="ERR2660591"
mkdir ${sample}
sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample="SRR16282278"
mkdir ${sample}
sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample="SRR11528297"
mkdir ${sample}
sbatch --job-name=$sample.ecc_caller_sra --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_sra.slurm

sample=ERR2660591
echo ${sample}
# count no lowq eccdnas
awk '$5 != "lowq"' ${sample}/${sample}.ecc_caller_out.details.txt | wc -l

# count split reads
wc -l ${sample}/${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}'

# count read numbers
samtools view -c -F 4 -F 2048 ${sample}/no_secondary.filtered.sorted.${sample}.bam

sample=SRR16282278
echo ${sample}
# count no lowq eccdnas
awk '$5 != "lowq"' ${sample}/${sample}.ecc_caller_out.details.txt | wc -l

# count split reads
wc -l ${sample}/${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}'

# count read numbers
samtools view -c -F 4 -F 2048 ${sample}/no_secondary.filtered.sorted.${sample}.bam

sample=SRR11528297
echo ${sample}
# count no lowq eccdnas
awk '$5 != "lowq"' ${sample}/${sample}.ecc_caller_out.details.txt | wc -l

# count split reads
wc -l ${sample}/${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}'

# count read numbers
samtools view -c -F 4 -F 2048 ${sample}/no_secondary.filtered.sorted.${sample}.bam