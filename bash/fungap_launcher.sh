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
    sbatch -p savio3 --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm
done < genomes_mapfile

while read genome; do
    echo ${genome}
    tail -2 ${genome}/fungap_out/logs/maker_ERR5875670_run1.log
done < genomes_mapfile