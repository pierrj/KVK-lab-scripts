#!/bin/bash

cd /global/scratch/users/pierrj/moryzae_virulence_project

while read genome; do
    echo $genome
    if [ -d "${genome}" ]; then
        rm -r ${genome}
    fi
    mkdir ${genome}
    cd ${genome}
        cp -r /global/scratch/users/pierrj/fungap_runs/gladieux_all/template_run/ .
    cd ..
done < strain_names

while read genome; do
    sbatch --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm
done < strain_names