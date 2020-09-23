#!/bin/bash
#SBATCH --job-name=word_cloud_wrapper
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/genes_on_eccs/neverever_wordcloud/
source activate pfamscan

SUBSET_GENE_IDS=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/genes_on_eccs/neverfoundoneccDNAs.geneids
ALL_GENE_IDS=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/genes_on_eccs/all_genenames.geneids
OUTPUT_NAME=G3_nevereverfound
PFAM_DIR=/global/scratch/users/pierrj/pfam
CDS_FASTA=/global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.CDS.fasta

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${SUBSET_GENE_IDS} -a ${ALL_GENE_IDS} -n ${OUTPUT_NAME} -p ${PFAM_DIR} -c ${CDS_FASTA}