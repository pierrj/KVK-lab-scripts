#!/bin/bash
#SBATCH --job-name=doesnt_matter
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

cd ${SLURM_SUBMIT_DIR}
cd ${sample}

## download reads from ENA, map reads and generate sorted bam files for ecc calling
/global/home/users/pierrj/git/bash/ena_generate_bam_file.sh -g ${genome_bwa} \
 -1 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_1.fastq.gz -2 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_2.fastq.gz \
 -s ${sample} -t ${SLURM_NTASKS} -m ${mapfile} -q ${READONE_LINK} -w ${READTWO_LINK}

## call putative eccDNA forming regions, using both uniquely mapped and multi-mapped reads
/global/home/users/pierrj/git/bash/ecc_caller_mapq0.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b uniq.filtered.sorted.${sample}.bam \
    -q multimapped.filtered.name_sorted.${sample}.bam

## assign confidence for eccDNA regions
/global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b no_secondary.filtered.sorted.${sample}.bam -r ${sample}.confirmedsplitreads.bed