#!/bin/bash

cd /global/scratch/users/pierrj/fungap_runs/wheat_blast/

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

while read genome; do
    echo ${genome}
    tail -2 ${genome}/fungap_out/logs/maker_ERR5875670_run1.log
done < genomes_mapfile


sbatch -p savio3 --ntasks-per-node=32 --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm\

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