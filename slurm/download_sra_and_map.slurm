#!/bin/bash
#SBATCH --job-name=download_sra_and_map.slurm
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd ${SLURM_SUBMIT_DIR}

module load bwa

if [ ! -f "${SAMPLE}.sra" ]; then
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${SLURM_NTASKS} -O . -t tmp ${SAMPLE}.sra
fi

READONE=${SAMPLE}.sra_1.fastq
READTWO=${SAMPLE}.sra_2.fastq

if [ ! -f "${SAMPLE}.bam" ]; then
    bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
    -t ${SLURM_NTASKS} ${REFERENCE} ${READONE} ${READTWO} -o ${SAMPLE}.preprocessed
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar SortSam \
        I=${SAMPLE}.preprocessed \
        O=${SAMPLE}.sorted \
        SORT_ORDER=coordinate
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar MarkDuplicates \
        I=${SAMPLE}.sorted \
        O=${SAMPLE}.bam \
        M=${SAMPLE}.marked_dup_metrics
    samtools index ${SAMPLE}.bam
fi

