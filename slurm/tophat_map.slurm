#!/bin/bash
#SBATCH --job-name=tophat_map
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
#cd /global/scratch/users/pierrj/references/
#bowtie2-build --threads 24 guy11_genome_baoetal2017.fasta guy11_genome_baoetal2017
module remove python/3.6
module load python/2.7
cd /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes
tophat -o 47_tophat -p 24 --b2-very-sensitive \
    --transcriptome-index=/global/scratch/users/pierrj/references/Guy11_augustus_raw_transcriptome --transcriptome-only \
    /global/scratch/users/pierrj/references/guy11_genome_baoetal2017 SRR298847.sra.fastq