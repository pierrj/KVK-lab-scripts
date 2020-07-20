#!/bin/bash
#SBATCH --job-name=get_srs_for_ltr
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/genome/2018_foster/

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames -s ERR2660591 -b sorted.mergedandpe.ERR2660591_bwamem.bam

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames -s G3_1A -b sorted.mergedandpe.G3_1A_bwamem.bam