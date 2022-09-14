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


genome=GCA_002924695.1_ASM292469v1

genome=GCA_905067075.2_PR003_contigs_polished

genome=GCA_905109835.1_Assembly_of_M.oryzae_isolate_KE017_genome

sbatch -p savio3 --ntasks-per-node=32 --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm


if [-f genomes_annotated_mapfile ]; then
    rm genomes_annotated_mapfile
fi

while read genome; do
    if [ -f ${genome}/fungap_out/fungap_out/fungap_out.gff3 ]; then
        echo ${genome} >> genomes_annotated_mapfile
    fi
done < genomes_mapfile


while read genome; do
    cp ${genome}/fungap_out/fungap_out/fungap_out_prot.faa /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/all_proteomes/${genome}_protein.fasta
done < genomes_annotated_mapfile



awk '{if ($3 > 90 && $3 != "Complete") {print substr($1, 0, length($1)-12)}}' genomes_for_busco_buscofied/batch_summary.txt > wheat_blast_busco_greater_than_90

squeue -u pierrj --format="%.100j" | tail -n +2 | awk '{print substr($1, 0,length($1)-11)}' > running_jobs

while read genome; do
    if [ ! -f ${genome}/fungap_out/fungap_out/fungap_out.gff3 ]; then
        if grep -Fxq "$genome" running_jobs
        then
            echo "$genome is already running"
        else
        sbatch -p savio3 --ntasks-per-node=32 --job-name=${genome}_run_fungap --export=genome=$genome /global/home/users/pierrj/git/slurm/run_fungap.slurm
        fi
    fi
done < wheat_blast_busco_greater_than_90


while read genome; do
    if [ -f ${genome}/fungap_out/fungap_out/fungap_out.gff3 ]; then
        echo ${genome} >> genomes_annotated_mapfile
    fi
done < wheat_blast_busco_greater_than_90

while read genome; do
    cp ${genome}/fungap_out/fungap_out/fungap_out_prot.faa /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/all_proteomes/${genome}_fungap_out_prot.faa
done < genomes_annotated_mapfile