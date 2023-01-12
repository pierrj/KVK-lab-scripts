#!/bin/bash
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

cd /global/scratch/users/pierrj/fungap_runs/wheat_blast/

# set up prelim files from templates for faster fungap runs
while read genome; do
    if [ -d "${genome}" ]; then
        rm -r ${genome}
    fi
    mkdir ${genome}
    cd ${genome}
        cp -r ../template_run/fungap_out/ .
    cd ..
done < genomes_mapfile

# submit fungap runs
while read genome; do
    sbatch --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm
done < genomes_mapfile

while read genome; do
    echo ${genome}
    tail -2 ${genome}/fungap_out/logs/maker_ERR5875670_run1.log
done < genomes_mapfile

# relaunch specific jobs
sbatch -p savio3 --ntasks-per-node=32 --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm

## to relaunch failed jobs due to busco download error

squeue -u pierrj --format="%.100j" | tail -n +2 | awk '{print substr($1, 0,length($1)-11)}' > running_jobs

while read genome; do
    if grep -Fxq "$genome" running_jobs
    then
        echo "$genome is already running"
    else
    sbatch -p savio3 --ntasks-per-node=32 --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm
    fi
done < genomes_mapfile