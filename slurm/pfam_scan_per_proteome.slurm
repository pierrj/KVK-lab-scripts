#!/bin/bash
#SBATCH --job-name=pfam_scan_per_proteome
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

DATA_OUTPUT_DIR=pfam_out

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

source activate /global/scratch/users/pierrj/conda_envs/pfam_scan

## run pfam_scan
pfam_scan.pl -cpu ${SLURM_NTASKS} -outfile pfamscan.out -dir ${PFAM_DIR} -fasta $PROTEIN_FASTA

## process and filter outputs
/global/scratch/users/pierrj/parse_pfam_domains/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p pfamscan.out -e 0.001 -o pfamscan.kparse.out

## generate table format
awk -v OFS='\t' '{print $1, "TRUE"}' pfamscan.kparse.out  > pfamscan.table