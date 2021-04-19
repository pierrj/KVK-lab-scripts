#!/bin/bash
#SBATCH --job-name=map_chpseq
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged

while read sample; do
    bwa mem -t ${SLURM_NTASKS} /global/scratch/users/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/guy11_genome_baoetal2017.fasta ${sample}.fastq -o ${sample}.sam
done < mapfile