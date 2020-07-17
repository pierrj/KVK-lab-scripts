#!/bin/bash
#SBATCH --job-name=ecc_dna_count_and_comparison
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

while read sample;
do
    cd ${sample}
    ecc_count=$(wc -l ecccaller_output.${sample}.bed | awk '{print $1}')
    read_count=$(samtools view -c guy11.repeats.sorted.mergedandpe.${sample}_bwamem.bam) ## should be number of reads mapped to scaffolds of interest
    genome_size=$(wc -l genomecoverage.mergedandpe.${sample}_bwamem.bed | awk '{print $1}') ## this is the new name for these, i think it is different for the old version of the pipelin
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}')
    echo -e $sample'\t'$total >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_count_comparison/ecc_count_comparison
    cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control

while read sample;
do
    cd ${sample}
    ecc_count=$(wc -l ecccaller_output.${sample}.bed | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam) ## should be number of reads mapped to scaffolds of interest
    genome_size=$(wc -l ${sample}.genomecoverage.bed | awk '{print $1}') ## this is the new name for these, i think it is different for the old version of the pipelin
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}')
    echo -e $sample'\t'$total >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_count_comparison/ecc_count_comparison
    cd ..
done < mapfile