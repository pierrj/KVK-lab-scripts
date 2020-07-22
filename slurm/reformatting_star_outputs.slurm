#!/bin/bash
#SBATCH --job-name=reformatting_STAR_outputs
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
for (( i = 47 ; i < 50; i++))
do
    cd /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes/${i}_star
    mv Aligned.out.sam ../${i}.Aligned.out.sam
    cd ..
    samtools view -b ${i}.Aligned.out.sam > ${i}.Aligned.out.bam
    samtools sort ${i}.Aligned.out.bam > sorted.${i}.Aligned.out.bam
    samtools index sorted.${i}.Aligned.out.bam
done