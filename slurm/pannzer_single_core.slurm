#!/bin/bash
#SBATCH --job-name=pannzer_single_core
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/scratch/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/scratch/users/pierrj/slurm_stderr/slurm-%j.out
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

## launch with
# proteome=something
# PROTEOMES_PATH=something

#  sbatch --job-name=$proteome.pannzer --mem=4000M -n 1 --export=ALL,proteome=$proteome,PROTEOMES_PATH=$PROTEOMES_PATH \
#     -p savio4_htc --qos=minium_htc4_normal --account=co_minium \
#     /global/home/users/pierrj/git/slurm/pannzer_single_core.slurm

source activate /global/scratch/users/pierrj/conda_envs/pannzer

cd $SLURM_SUBMIT_DIR

protein_count=$(grep gene GO.${proteome}.out | awk '{print $1}' | sort | uniq | wc -l)
b=7000

## check to see if run is already done, sometimes they dont fully finish
if (( protein_count < b )); then
    mv GO.${proteome}.out GO.${proteome}.out.old 
    ## run pannzer off of the pannzer server
    /global/scratch/users/pierrj/conda_envs/pannzer/bin/python \
        /global/scratch/users/pierrj/pannzer2/SANSPANZ.3/runsanspanz.py \
        -R -o ",DE.${proteome}.out,GO.${proteome}.out,anno.${proteome}.out" -s "Pyricularia oryzae" < ${PROTEOMES_PATH}/${proteome}
fi