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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper/anyoverlap_normalization

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
SAMPLE=G3
SRA_LIST=guy11_sra_accessions
THREADS=${SLURM_NTASKS}
LIB_TYPE=2
GFF_FILE=/global/scratch/users/pierrj/references/guy11_fungap.gff3
ECCDNA_MAPFILE=G3.eccmapfile
SAMPLE_MAPFILE=G3.sample_mapfile

module load STAR

/global/home/users/pierrj/git/bash/ecc_dna_vs_expression.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n a

/global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper/multiplicative_normalization

/global/home/users/pierrj/git/bash/ecc_dna_vs_expression.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n g