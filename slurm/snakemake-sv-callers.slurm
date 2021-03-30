#!/bin/bash
#SBATCH --job-name=snakemake-sv-callers
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


REFERENCE="moryzae_70-15_ref_with_mito.fasta"
SAMPLE=ERR4757126
GENOME_DB=/global/scratch/users/pierrj/sv-callers_runs/guy11_v_AG038/data/fasta/moryzae_70-15_ref_with_mito.fasta
READONE=${SAMPLE}.sra_1.fastq
READTWO=${SAMPLE}.sra_2.fastq
SNAKEMAKE_DIR=/global/scratch/users/pierrj/sv-callers_runs/guy11_v_AG038
GIT_REPO=/global/home/users/pierrj/git
ANALYSIS_YAML=snakemake/analysis.yaml
SAMPLES_CSV=snakemake/samples.csv

cd ${SNAKEMAKE_DIR}/data/fasta

if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx ${REFERENCE}
fi

if [ ! -f "${REFERENCE}.bwt" ]; then
    bwa index ${REFERENCE}
fi

cd ${SNAKEMAKE_DIR}/data/bam/${SAMPLE}

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

source ~/.bashrc
conda activate wf

module purge

cd ${GIT_REPO}

git pull

cd ${SNAKEMAKE_DIR}

cp ${GIT_REPO}/${ANALYSIS_YAML} analysis.yaml
cp ${GIT_REPO}/${SAMPLES_CSV} samples.csv

snakemake --use-conda --latency-wait 30 --jobs 14 --cluster "sbatch --partition=savio2 --account=fc_kvkallow --qos=savio_normal --nodes=1 --cpus-per-task={threads} --time=72:00:00 --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out --job-name=smk.{rule} --mail-user=pierrj@berkeley.edu --mail-type=ALL"