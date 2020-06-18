#!/bin/bash
#SBATCH --job-name=anygenome_coverageplot
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/anygenome_coverageplot_guy11_test

# /global/home/users/pierrj/git/bash/generate_bam_file_genomecov.sh -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_tweaked_bwa \
#    -1 /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R1.fastq -2 /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R2.fastq \
#    -s CT_1A -t 20

# /global/home/users/pierrj/git/bash/generate_coverageplot.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames -s CT_1A -c CT_1A.genomecoverage.bed

/global/home/users/pierrj/git/bash/call_ecc_regions_flags_test.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames\
    -s CT_1A -t 20 -c CT_1A.genomecoverage.bed -b /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/sorted.mergedandpe.CT_1A_bwamem.bam