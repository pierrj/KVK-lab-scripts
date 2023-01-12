#!/bin/bash
#SBATCH --job-name=download_sra_and_map.slurm
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
#Copyright (c) 2023 Pierre Michel Joubert
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

module load bwa

## download data from sra accession, has to be paired end
if [ ! -f "${SAMPLE}.sra" ]; then
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${SLURM_NTASKS} -O . -t tmp ${SAMPLE}.sra
fi

READONE=${SAMPLE}.sra_1.fastq
READTWO=${SAMPLE}.sra_2.fastq

## map, sort, mark duplicates
if [ ! -f "${SAMPLE}.bam" ]; then
    bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
    -t ${SLURM_NTASKS} ${REFERENCE} ${READONE} ${READTWO} -o ${SAMPLE}.preprocessed
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar SortSam \
        I=${SAMPLE}.preprocessed \
        O=${SAMPLE}.sorted \
        SORT_ORDER=coordinate
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar MarkDuplicates \
        I=${SAMPLE}.sorted \
        O=${SAMPLE}.bam \
        M=${SAMPLE}.marked_dup_metrics
    samtools index ${SAMPLE}.bam
fi

