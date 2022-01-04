#!/bin/bash
#SBATCH --job-name=moryzae_sv_calling_manta_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/data

REFERENCE=moryzae_70-15_ref_with_mito.fasta
SAMPLE=ERR4757126
GENOME_DB=moryzae_70-15_ref_with_mito.fasta
READONE=${SAMPLE}.sra_1.fastq
READTWO=${SAMPLE}.sra_2.fastq

if [ ! -f "${SAMPLE}.sra_1.fastq" ]; then
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${SLURM_NTASKS} -O . -t tmp ${SAMPLE}.sra
fi

if [ ! -f "${REFERENCE}.amb" ]; then
    bwa index ${REFERENCE}
fi

if [ ! -f "${SAMPLE}.bam" ]; then
    bwa mem -t ${SLURM_NTASKS} ${GENOME_DB} ${READONE} ${READTWO} -o ${SAMPLE}.preprocessed
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar SortSam \
        I=${SAMPLE}.preprocessed \
        O=${SAMPLE}.sorted \
        SORT_ORDER=coordinate
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar MarkDuplicates \
        I=${SAMPLE}.sorted \
        O=${SAMPLE}.bam \
        M=${SAMPLE}.marked_dup_metrics
fi

if [ ! -f "${SAMPLE}.bam.bai" ]; then
    samtools index ${SAMPLE}.bam
fi

cd /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/delly

module purge

source activate /global/scratch/users/pierrj/conda_envs/manta

/global/scratch/users/pierrj/manta-1.6.0.centos6_x86_64/bin/configManta.py \
    --runDir manta_out \
    --reference ../data/${REFERENCE} \
    --bam ../${SAMPLE}.bam

/manta_out/runWorkflow.py \
    --quiet \
    -m local \
    -j ${SLURM_NTASKS}

## no filtering