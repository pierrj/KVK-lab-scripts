#!/bin/bash
#SBATCH --job-name=genomecov_100kbbins
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
cd /global/scratch/users/pierrj/references
bedtools intersect -a guy11_genome_baoetal2017.100kbbins -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/sorted.mergedandpe.G3_1A_bwamem.bam -c -sorted -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > /global/scratch/users/pierrj/eccDNA/pipeline_tests/codingregions_comparison/ecc_rawreads_perk100kb