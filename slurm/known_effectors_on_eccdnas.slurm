#!/bin/bash
#SBATCH --job-name=known_effectors_on_eccdnas
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

# effector location
effector_name=AvrPita1
effector_chrom=MQOP01000011.1
effector_start=675417
effector_end=677485
effector_gene_name=gene_10084

# see if any eccdnas fully encompass it
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

# found?
## yes

# effector location
effector_name=AvrPita3
effector_chrom=MQOP01000016.1
effector_start=368546
effector_end=370280
effector_gene_name=gene_11229

# see if any eccdnas fully encompass it
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

# found?
## yes


# effector location
effector_name=AvrPi9
effector_chrom=MQOP01000008.1
effector_start=2277217
effector_end=2277733
effector_gene_name=gene_09052

# see if any eccdnas fully encompass it
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

# found?
## yes

# effector location
effector_name=AvrPi54
effector_chrom=MQOP01000009.1
effector_start=815984
effector_end=816445
effector_gene_name=gene_09359

# see if any eccdnas fully encompass it
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

# found?
## yes

# effector location
effector_name=AvrPiz-t
effector_chrom=MQOP01000012.1
effector_start=880675
effector_end=881001
effector_gene_name=gene_10438

# see if any eccdnas fully encompass it
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

# found?
## yes

# effector location
effector_name=Pwl4
effector_chrom=MQOP01000002.1
effector_start=4264949
effector_end=4265365
effector_gene_name=gene_02871

# see if any eccdnas fully encompass it
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

# found?
## yes