#!/bin/bash
#SBATCH --job-name=gene_portions_moryzae
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/script_test/This_Study=Rice_blast=Mycelium=WT=None


ORGANISM_NAME="This_Study=Rice_blast=Mycelium=WT=None"
GENE_GFF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
SAMPLE_MAPFILE=mapfile
ECC_DIR=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/
GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
COPIA_FILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_copia_locs.bed
GYPSY_FILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_gypsy_locs.bed
CONTIGNAMES_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames


/global/home/users/pierrj/git/bash/gene_portions_enrichment.sh -o ${ORGANISM_NAME} \
    -g ${GENE_GFF} -s ${SAMPLE_MAPFILE} -d ${ECC_DIR} \
    -f ${GENOME_FILE} -c ${COPIA_FILE} -y ${GYPSY_FILE} -n ${CONTIGNAMES_FILE}

mv observed_averages_large_all_samples moryzae.observed_averages_large_all_samples
mv observed_averages_micro_all_samples moryzae.observed_averages_micro_all_samples