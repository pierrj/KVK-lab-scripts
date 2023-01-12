#!/bin/bash
#SBATCH --job-name=fasttree_scos
#SBATCH --partition=savio3_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=48:00:00
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


PROJECT_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all
DATE=Sep08

cd $PROJECT_DIR

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

SCO_DIR=$PROJECT_DIR/orthofinder/Results_${DATE}/WorkingDirectory/OrthoFinder/Results_${DATE}/Single_Copy_Orthologue_Sequences
basename=$(basename ${SCO_DIR})
MSA=${basename}_msas_combined.trimmed.fasta

# make phylogeny from msa
fasttree -gamma -out ${basename}_msas_combined.trimmed.fasttree ${basename}_msas_combined.trimmed.fasta