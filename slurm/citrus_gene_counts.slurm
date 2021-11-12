#!/bin/bash
#SBATCH --job-name=citrus_gene_counts
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/citrus_rnaseq

THREADS=24
GENOME_FASTA=GCF_000493195.1_Citrus_clementina_v1.0_genomic.fna
GFF_FILE=GCF_000493195.1_Citrus_clementina_v1.0_genomic.gff
SRA_LIST=c_bergamia_runs

# get genome and unzip
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/493/195/GCF_000493195.1_Citrus_clementina_v1.0/GCF_000493195.1_Citrus_clementina_v1.0_genomic.fna.gz
gunzip GCF_000493195.1_Citrus_clementina_v1.0_genomic.fna.gz

# get gff and unzip
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/493/195/GCF_000493195.1_Citrus_clementina_v1.0/GCF_000493195.1_Citrus_clementina_v1.0_genomic.gff.gz
gunzip GCF_000493195.1_Citrus_clementina_v1.0_genomic.gff.gz

if [ -f "${SRA_LIST}" ]; then
    rm ${SRA_LIST}
fi

## add run ids
echo ERR760726 >> ${SRA_LIST} ## aka SAMEA2612286
echo ERR760727 >> ${SRA_LIST} ## aka SAMEA2612285

genome_fasta_basename=$(basename ${GENOME_FASTA})

if [ -d "${genome_fasta_basename}_starindex" ]; then
    rm -r ${genome_fasta_basename}_starindex
fi

mkdir ${genome_fasta_basename}_starindex

module load STAR

## index genome for STAR
STAR --runThreadN ${THREADS} --runMode genomeGenerate --genomeDir ${genome_fasta_basename}_starindex \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${GFF_FILE} \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11 \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent


## map
while read SRA; do
        /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
        /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
        STAR --runThreadN ${THREADS} \
            --genomeDir ${genome_fasta_basename}_starindex \
            --readFilesIn ${SRA}.sra_1.fastq ${SRA}.sra_2.fastq \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${SRA}. \
            --quantMode GeneCounts
done < ${SRA_LIST}

