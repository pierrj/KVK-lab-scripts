#!/bin/bash
#SBATCH --job-name=eccdna_circrna
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression


## trimmomatic or something maybe?

SRA=SRR8842990

bwa mem –T 19 -t ${SLURM_NTASKS} guy11_genome_baoetal2017.fasta ${SRA}.sra_1.fastq ${SRA}.sra_2.fastq > ${SRA}.sam

perl /global/scratch/users/pierrj/ciri2/CIRI_v2.0.5/CIRI_v2.0.5.pl -T ${SLURM_NTASKS} -I ${SRA}.sam -O ${SRA}.ciri_out -F guy11_genome_baoetal2017.fasta