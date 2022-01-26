#!/bin/bash

cd /global/scratch/users/pierrj/fungap_runs/gladieux_all/

while read genome; do
    if [ -d "${genome}" ]; then
        rm -r ${genome}
    fi
    mkdir ${genome}
    cd ${genome}
        cp -r ../template_run/fungap_out/ .
    cd ..
done < genomes_mapfile

while read genome; do
    sbatch --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm
done < genomes_mapfile