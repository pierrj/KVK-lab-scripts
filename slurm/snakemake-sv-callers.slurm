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


REFERENCE="guy11_genome_baoetal2017.fasta"
SAMPLE=ERR4757126
GENOME_DB=/global/scratch/users/pierrj/sv-callers_runs/guy11_v_AG038/data/fasta/guy11_genome_baoetal2017.fasta
READONE=${SAMPLE}.sra_1.fastq
READTWO=${SAMPLE}.sra_2.fastq
SNAKEMAKE_DIR=/global/scratch/users/pierrj/sv-callers_runs/guy11_v_AG038
GIT_REPO=/global/home/users/pierrj/git
ANALYSIS_YAML=snakemake/analysis.yaml

cd ${SNAKEMAKE_DIR}/data/fasta

if [ ! -f "${REFERENCE}.fai" ]; then
    samtools index ${REFERENCE}
fi

if [ ! -f "${REFERENCE}.bwt" ]; then
    bwa index ${REFERENCE}
fi

cd ${SNAKEMAKE_DIR}/data/bam/${SAMPLE}

if [ ! -f "${SAMPLE}_bwamem.bam" ]; then
    bwa mem -t ${SLURM_NTASKS} ${GENOME_DB} ${READONE} ${READTWO} -o ${SAMPLE}_bwamem.preprocessed
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar SortSam \
        I=${SAMPLE}_bwamem.preprocessed \
        O=${SAMPLE}_bwamem.sorted \
        SORT_ORDER=coordinate
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar MarkDuplicates \
        I=${SAMPLE}_bwamem.sorted \
        O=${SAMPLE}_bwamem.bam \
        M=${SAMPLE}_bwamem.marked_dup_metrics
fi

source ~/.bashrc
conda activate wf

cd ${GIT_REPO}

git pull

cd ${SNAKEMAKE_DIR}

cp ${GIT_REPO}/${ANALYSIS_YAML} analysis.yaml

snakemake --use-conda --latency-wait 30 --jobs 14 --cluster "sbatch --partition=savio2 --account=fc_kvkallow --qos=savio_normal --nodes=1 --cpus-per-task={threads} --time=72:00:00 --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out --job-name=smk.{rule} --mail-user=pierrj@berkeley.edu --mail-type=ALL"