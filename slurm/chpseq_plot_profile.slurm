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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged

source activate deeptools

# while read biorep treatment; do
#     samtools sort ${treatment}_${biorep}.sam.bam > sorted.${treatment}_${biorep}.sam.bam
#     samtools index sorted.${treatment}_${biorep}.sam.bam
# done < biorep_treament

# while read biorep treatment; do
#     samtools sort ${treatment}_${biorep}.sam.bam > sorted.${treatment}_${biorep}.sam.bam
#     samtools index sorted.${treatment}_${biorep}.sam.bam
# done < control_biorep_treament

cd ../deeptools

while read biorep treatment; do
    bamCompare -b1 ../samples_merged/sorted.${treatment}_${biorep}.sam.bam -b2 ../samples_merged/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
    computeMatrix scale-regions -S log2ratio_${treatment}_${biorep}.bw \
                              -R ../sorted.uniq.all.confirmedsplitreads.bed ../shuffled_eccs.bed \
                              --beforeRegionStartLength 300 \
                              --regionBodyLength 1000 \
                              --afterRegionStartLength 300 \
                              --skipZeros -o /tmp/${treatment}_${biorep}.mat.gz 
    plotProfile -m /tmp/${treatment}_${biorep}.mat.gz \
              -out /tmp/${treatment}_${biorep}.png \
              --numPlotsPerRow 1 \
              --plotTitle "${treatment} ${biorep} eccDNAs"
done < biorep_treament

mv /tmp/* /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/