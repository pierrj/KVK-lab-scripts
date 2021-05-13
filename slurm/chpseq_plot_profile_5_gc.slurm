#!/bin/bash
#SBATCH --job-name=chpseq_plot_profile
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/tracks5

source activate deeptools

# while read biorep treatment; do
#     samtools sort ${treatment}_${biorep}.sam.bam > sorted.${treatment}_${biorep}.sam.bam
#     samtools index sorted.${treatment}_${biorep}.sam.bam
# done < biorep_treament

# while read biorep treatment; do
#     samtools sort ${treatment}_${biorep}.sam.bam > sorted.${treatment}_${biorep}.sam.bam
#     samtools index sorted.${treatment}_${biorep}.sam.bam
# done < control_biorep_treament

biorep=1

computeMatrix scale-regions -S ../gc/guy11_genome_baoetal2017.50bpwindows.bw \
                            -R genomic_baseline all_genes common_genes neverfound_genes \
                            --beforeRegionStartLength 300 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 300 \
                            --skipZeros -o gc.mat.gz

plotProfile -m gc.mat.gz \
            -out gc.png \
            --numPlotsPerRow 1 \
            --plotTitle "GC eccDNAs"