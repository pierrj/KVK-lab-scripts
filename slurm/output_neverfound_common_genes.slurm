#!/bin/bash
#SBATCH --job-name=output_neverfound_common_genes
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/compare_bioreps

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav/G3.all_eccdnas_mapfile
OUTPUT_NAME=G3
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile

# output never (eccdna-absent) and common (eccdna-associated) genes using encompassing split reads
/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -n ${OUTPUT_NAME} \
    -b ${GENE_BEDFILE} -i ${SAMPLE_MAPFILE} -e ${ECCDNA_MAPFILE}