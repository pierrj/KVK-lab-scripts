#!/bin/bash
#SBATCH --job-name=percent_ltr_eccs_circularome
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

if [ -f "ltr_eccs_count" ]; then
    rm ltr_eccs_count
fi

if [ -f "ltr_eccs_splitreads_count" ]; then
    rm ltr_eccs_splitreads_count
fi

if [ -f "ecc_region_count" ]; then
    rm ecc_region_count
fi

if [ -f "junction_sr_count" ]; then
    rm junction_sr_count
fi

# count number of ltr eccs, ltr junction split reads, ecc regions and total junction split reads
while read sample; do
cd ${sample}
    wc -l ltr_eccs | awk '{print $1}' >> ../ltr_eccs_count
    wc -l ltr_eccs_splitreads | awk '{print $1}' >> ../ltr_eccs_splitreads_count
    wc -l ${sample}.ecc_caller_out.details.nolowq.txt | awk '{print $1}' >> ../ecc_region_count
    wc -l ${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}' >> ../junction_sr_count
cd ..
done < mapfile

# print counts
awk '{sum+=$1} END {print sum}' ltr_eccs_count
awk '{sum+=$1} END {print sum}' ltr_eccs_splitreads_count
awk '{sum+=$1} END {print sum}' ecc_region_count
awk '{sum+=$1} END {print sum}' junction_sr_count
