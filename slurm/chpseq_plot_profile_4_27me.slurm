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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/tracks4

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

treatment=H3K27me3
bamCompare -b1 ../samples_merged/sorted.${treatment}_${biorep}.sam.bam -b2 ../samples_merged/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
computeMatrix scale-regions -S log2ratio_${treatment}_${biorep}.bw \
                            -R  genomic_baseline all_genes confident_ltr_tes ltr_eccdnas micro_dnas large_eccdnas methylated_genes nonmethylated_genes \
                            --beforeRegionStartLength 300 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 300 \
                            --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.png \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs"