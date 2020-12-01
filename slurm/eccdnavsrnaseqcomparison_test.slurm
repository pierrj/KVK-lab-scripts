#!/bin/bash
#SBATCH --job-name=eccdnavsrnaseqcomparison_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper

module load STAR

/global/home/users/pierrj/git/bash/ecc_dna_vs_expression.sh -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta \
    -s G3 -a guy11_sra_accessions -t ${SLURM_NTASKS} -l 2 -f /global/scratch/users/pierrj/references/guy11_fungap.gff3 \
    -e G3.eccmapfile -m G3.sample_mapfile