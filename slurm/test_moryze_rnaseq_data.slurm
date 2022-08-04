#!/bin/bash
#SBATCH --job-name=test_moryze_rnaseq_data
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/kamoun_lab_rnaseq

THREADS=20
GENOME_FASTA=guy11_plus_orsa_w_organelles.fasta
GFF_FILE=guy11_plus_orsa.gff

SRA_LIST=italian_samples

if [ -f "${SRA_LIST}" ]; then
    rm ${SRA_LIST}
fi

## add run ids
echo ERR2593104 >> ${SRA_LIST}
echo ERR2593112 >> ${SRA_LIST}
echo ERR2593113 >> ${SRA_LIST}
echo ERR2593117 >> ${SRA_LIST}
echo ERR2593119 >> ${SRA_LIST}
echo ERR2593106 >> ${SRA_LIST}
echo ERR2593108 >> ${SRA_LIST}
echo ERR2593110 >> ${SRA_LIST}

genome_fasta_basename=$(basename ${GENOME_FASTA})

if [ -d "${genome_fasta_basename}_starindex" ]; then
    rm -r ${genome_fasta_basename}_starindex
fi

mkdir ${genome_fasta_basename}_starindex

module load STAR

## index genome for STAR
## maybe not GFF next time?
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