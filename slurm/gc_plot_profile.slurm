#!/bin/bash
#SBATCH --job-name=gc_plot_profile
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/gc

source activate deeptools

computeMatrix scale-regions -S log2ratio_gc.bw \
                            -R ../sorted.uniq.all.confirmedsplitreads.bed ../shuffled_eccs.bed \
                            --beforeRegionStartLength 300 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 300 \
                            --skipZeros -o gc.mat.gz

plotProfile -m gc.mat.gz \
            -out gc.png \
            --numPlotsPerRow 1 \
            --plotTitle "GC eccDNAs"