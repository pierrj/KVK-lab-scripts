#!/bin/bash
#SBATCH --job-name=ecc_finder
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
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

## LAUNCH JOB WITH
# while read sample; do 
#     sbatch --job-name=$sample.ecc_finder --export=sample=$sample /global/home/users/pierrj/git/slurm/ecc_finder.slurm
# done < /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/mapfile


sample=G3_3A
genome=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
reads_1=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}_R1.fastq
reads_2=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}_R2.fastq

cd /global/scratch/users/pierrj/ecc_finder

source activate /global/scratch/users/pierrj/conda_envs/ecc_finder_mine

export PATH=/global/scratch/users/pierrj/conda_envs/ecc_finder_mine/bin:$PATH

python ecc_finder.py map-sr ${genome} \
            ${reads_1} \
            ${reads_2} \
            -r ${genome} \
            -t ${SLURM_NTASKS} -o ecc_finder_${sample}