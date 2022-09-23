#!/bin/bash
#SBATCH --job-name=ecc_finder
#SBATCH --partition=savio2
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


sample=SRR6315399

DIR=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/ecc_finder_${sample}

if [ -d "$DIR" ]; then
    rm -r $DIR
fi

mkdir $DIR

cd $DIR

mkdir align_files

cd align_files

genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"

source activate /global/scratch/users/pierrj/conda_envs/ecc_finder_pre_steps

module purge

cutadapt -j ${SLURM_NTASKS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --nextseq-trim=20 -o ecc.sr.R1.fastq.gz \
    -p ecc.sr.R2.fastq.gz \
    /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/${sample}.sra_1.fastq \
    /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/${sample}.sra_2.fastq


# fastp --thread ${SLURM_NTASKS} -i /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/${sample}.sra_1.fastq \
#             -I /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/${sample}.sra_2.fastq \
#             -o ecc.sr.R1.fastq.gz \
#             -O ecc.sr.R2.fastq.gz \
#             -h ecc.sr.html

bwa mem -t ${SLURM_NTASKS} ${genome_bwa} \
            ecc.sr.R1.fastq.gz \
            ecc.sr.R2.fastq.gz > ecc.sr.sam 2> ecc.sr.sam.log

samtools sort -@ ${SLURM_NTASKS} -n -o ecc.sr.bam ecc.sr.sam

rm ecc.sr.sam

touch ecc.sr.sam

## transfer to Mac pro then run ecc_finder there since it doesn't work on cluster
# python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 -l 0 --min-read 0 -p 1.0 -o test_human >& run.test_human.${params}.out