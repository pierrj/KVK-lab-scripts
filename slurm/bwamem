#!/bin/bash
#SBATCH --job-name=bwamem.singleend
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/G3_1A_tests
bwa mem -t 24 guy11_genome_baoetal2017_bwa forsingleendmapping.G3_1A.fastq -o singleend.G3_1A_bwamem.sam
