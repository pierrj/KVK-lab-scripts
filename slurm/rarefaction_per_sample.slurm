#!/bin/bash
#SBATCH --job-name=rarefaction_per_sample
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
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/rarefaction_tech_reps/${SAMPLE}

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"

for i in $(seq 0.1 0.1 1.0); do

    sample=${i}

    if [ -d "${i}" ]; then
        rm -r ${i}
    fi

    mkdir ${i}

    cd ${i}

    samtools view -s ${i} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/uniq.filtered.sorted.${SAMPLE}.bam > uniq.filtered.sorted.${sample}.bam
    samtools view -s ${i} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/multimapped.filtered.name_sorted.${SAMPLE}.bam > multimapped.filtered.name_sorted.${sample}.bam
    samtools view -s ${i} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/no_secondary.filtered.sorted.${SAMPLE}.bam > no_secondary.filtered.sorted.${sample}.bam

    ## call putative eccDNA forming regions, using both uniquely mapped and multi-mapped reads
    /global/home/users/pierrj/git/bash/ecc_caller_mapq0.sh -m ${mapfile} \
        -s ${sample} \
        -t ${SLURM_NTASKS} \
        -b uniq.filtered.sorted.${sample}.bam \
        -q multimapped.filtered.name_sorted.${sample}.bam

    ## assign confidence for eccDNA regions
    /global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
        -s ${sample} \
        -t ${SLURM_NTASKS} \
        -b no_secondary.filtered.sorted.${sample}.bam \
        -r ${sample}.confirmedsplitreads.bed

    cd ..
done 