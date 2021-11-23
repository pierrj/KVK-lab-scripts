#!/bin/bash
#SBATCH --job-name=junction_split_reads_per_sample
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


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

# generate table
output_table=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/srs_per_reads

if [ -f "${output_table}" ]; then
    rm ${output_table}
fi

while read sample; do
    cd ${sample}
        # count reads
        read_count=$(samtools view -c -F 4 -F 2048 no_secondary.filtered.sorted.${sample}.bam)
        # count srs
        sr_count=$(wc -l ${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}')
        echo -e ${sample}'\t'${read_count}'\t'${sr_count} >> /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/srs_per_reads
    cd ..
done < mapfile


awk '{print $3/$2*1000000}' ${output_table} | awk '{SUM+=$1} END {print SUM/NR}'