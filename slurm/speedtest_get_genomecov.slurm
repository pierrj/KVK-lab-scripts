#!/bin/bash
#SBATCH --job-name=speedtest_get_genomecov
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/speedtest_get_genomecov

sample=RC_1A
bamfile=/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A.sorted.mergedandpe.bwamem.bam

date +"%T"

source activate mosdepth

mosdepth ${sample} ${bamfile} -t 3 -x -b 1

date +"%T"

samtools depth -a ${bamfile} > ${sample}.samtoolsdepth.bed

date +"%T"

bedtools genomecov -d -ibam ${bamfile} > ${sample}.genomecoverage.bed

date +"%T"