#!/bin/bash
#SBATCH --job-name=repeatmasker_rice_blast
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

PROJECT_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap

if [ ! -d "$PROJECT_DIR/repeatmasked" ] 
then
    mkdir $PROJECT_DIR/repeatmasked
fi

cd $PROJECT_DIR/repeatmasked

RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
GENOMES_PATH=/global/scratch/users/pierrj/fungap_runs/gladieux_all/genomes_to_annotate
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/genomes_mapfile

module load repeatmasker/4.1.1

## repeatmask all genomes
while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}_genomic.fna
done < ${MAPFILE}

while read genome; do
    awk '{print $11}' ${genome}_genomic.fna.out | paste ${genome}_genomic.fna.out.gff - > ${genome}.w_class.gff
done < ${MAPFILE}