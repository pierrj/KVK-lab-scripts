#!/bin/bash
#SBATCH --job-name=test_map_moryzae
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/rnaseq/moryzae/2021_oses-ruiz

module load STAR

THREADS=24
GENOME_FASTA=guy11_genome_baoetal2017.fasta
GFF_FILE=guy11_fungap_out_12_28_20.gff3

genome_fasta_basename=$(basename ${GENOME_FASTA})

# mkdir ${genome_fasta_basename}_starindex

# ## index genome for STAR
# STAR --runThreadN ${THREADS} --runMode genomeGenerate --genomeDir ${genome_fasta_basename}_starindex \
#     --genomeFastaFiles ${GENOME_FASTA} \
#     --sjdbGTFfile ${GFF_FILE} \
#     --sjdbOverhang 100 \
#     --genomeSAindexNbases 11 \
#     --sjdbGTFtagExonParentTranscript ID \
#     --sjdbGTFtagExonParentGene Parent

SRA=ERR5875672

/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
STAR --runThreadN ${THREADS} \
    --genomeDir ${genome_fasta_basename}_starindex \
    --readFilesIn ${SRA}.sra_1.fastq ${SRA}.sra_2.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${SRA}. \
    --quantMode GeneCounts