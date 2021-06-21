#!/bin/bash
#SBATCH --job-name=circos_plot
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot


bedtools makewindows -g guy11_genome_baoetal2017.chromsizes -w 1000 > guy11_genome_baoetal2017.1000bpwindows

bedtools coverage -a guy11_genome_baoetal2017.1000bpwindows -b G3_1A.ecc_caller_out.splitreads.bed -sorted -g guy11_genome_baoetal2017.chromsizes | awk -v OFS='\t' '{print $1, $2, $3, $4}' > eccdna_genomecoverage_G3_1A

## te bed track to lines

## gene bed track to lines

## some sort of density plot for eccdnas/splitreads

## https://jokergoo.github.io/circlize_book/book/introduction.html

## http://felixfan.github.io/circos/