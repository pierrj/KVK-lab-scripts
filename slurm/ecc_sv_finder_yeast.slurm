#!/bin/bash
#SBATCH --job-name=ecc_sv_finder_yeast
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/structural_variation/yeast

module load gnuplot
module load imagemagick

## copy yeast genome fastas

cp /global/scratch/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script/yeast_test/EC1118_yeast.fasta .
cp /global/scratch/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script/yeast_test/M22_yeast.fasta .


## get ecc svs for EC1118 as reference and M22 as query
/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r EC1118_yeast.fasta -q M22_yeast.fasta -t ${SLURM_NTASKS} -n ec1118_v_m22 \
    -p 1 -o ec1118_v_m22_out -m ec1118_v_m22_tmp

## get ecc svs for M22 as reference and EC1118 as query
/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r M22_yeast.fasta -q EC1118_yeast.fasta -t ${SLURM_NTASKS} -n m22_v_ec1118 \
    -p 1 -o m22_v_ec1118_out -m m22_v_ec1118_tmp