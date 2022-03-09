#!/bin/bash
#SBATCH --job-name=coverage_for_pca
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/qc/pca_coverage_tech_reps

CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes


bedtools makewindows -g ${CHROM_SIZES} -w 1000 > guy11_genome_baoetal2017.1000bpwindows


while read sample
do
# junction split reads corresponding to each window
bedtools coverage -a guy11_genome_baoetal2017.1000bpwindows \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.genomecoverage.bed
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile