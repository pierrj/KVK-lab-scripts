#!/bin/bash
#SBATCH --job-name=id_known_effectors
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

cd /global/scratch/users/pierrj/diamond/guy11_avr_locs

# proteome from fungap
../diamond makedb --in guy11_fungap_out_12_28_20_prot.faa -d guy11_fungap_12_28_20_prot

# avrs from Kim et al 2019
../diamond blastp -d guy11_fungap_12_28_20_prot -q M.oryzae.AVRs.fasta -o matches.tsv


## filter out identity less than 50%, then make sure each fungap gene only has the best match assigned to it
awk '$3>=50.0' matches.tsv | awk '$2 != "gene_09984.t1" && $1 != "AvrPita2" && $1 != "Pwl3" && !($1 == "AvrPita3" && $2 == "gene_10084.t1")' \
 | awk -v OFS='\t' '{print $1, $2}' > known_effectors_fungap_names

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common

cp /global/scratch/users/pierrj/diamond/guy11_avr_locs/known_effectors_fungap_names .