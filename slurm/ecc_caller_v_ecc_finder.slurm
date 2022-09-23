#!/bin/bash
#SBATCH --job-name=ecc_caller_v_ecc_finder
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human_ecc_finder

sample_mapfile=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/mapfile
mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"

sample=SRR6315399

cd ${sample}

/global/home/users/pierrj/git/bash/ecc_caller_mapq0_nomatchlength.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/uniq.filtered.sorted.${sample}.bam \
    -q /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/multimapped.filtered.name_sorted.${sample}.bam

# call confidence
/global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/no_secondary.filtered.sorted.${sample}.bam \
    -r ${sample}.unique.confirmedsplitreads.bed

mv ${sample}.ecc_caller_out.details.txt ${sample}.ecc_caller_out.no_length_filter.details.txt