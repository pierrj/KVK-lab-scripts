#!/bin/bash
#SBATCH --job-name=rarefaction_per_tech_rep
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/rarefaction_tech_reps/

MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

if [ -f "rarefaction_table" ]; then
    rm rarefaction_table
fi

while read sample; do
    if [ -f "${sample}_column" ]; then
        rm ${sample}_column
    fi

    for i in $(seq 0.1 0.1 1.0); do
        echo ${sample} >> ${sample}_column
    done

    if [ -f "${sample}_rarefaction_percent" ]; then
        rm ${sample}_rarefaction_percent
    fi

    for i in $(seq 0.1 0.1 1.0); do
        echo ${i} >> ${sample}_rarefaction_percent
    done

    if [ -f "${sample}_read_counts" ]; then
        rm ${sample}_read_counts
    fi

    for i in $(seq 0.1 0.1 1.0); do
        samtools view -c -F 4 -F 2048 ${sample}/${i}/no_secondary.filtered.sorted.${i}.bam >> ${sample}_read_counts
    done

    if [ -f "${sample}_ecc_counts" ]; then
        rm ${sample}_ecc_counts
    fi

    for i in $(seq 0.1 0.1 1.0); do
        sort ${sample}/${i}/${i}.ecc_caller_out.splitreads.bed | uniq | wc -l >> ${sample}_ecc_counts
    done

    paste ${sample}_column ${sample}_rarefaction_percent ${sample}_read_counts ${sample}_ecc_counts > ${sample}_rarefaction_curve

    cat ${sample}_rarefaction_curve >> rarefaction_table

done < $MAPFILE