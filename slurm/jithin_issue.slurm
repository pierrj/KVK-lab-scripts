#!/bin/bash
#SBATCH --job-name=jithin_issue
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/jithin_issue/SRR6315431

SAMPLE=SRR6315431
THREADS=${SLURM_NTASKS}

/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SAMPLE}.sra

export ECC_CALLER_PYTHON_SCRIPTS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/jithin_issue/git/ecc_caller/python_scripts/
export ECC_CALLER_PICARD=/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

bwa index -p hg19_bwa hg19.fa

GENOME_BWA=hg19_bwa

grep '>' hg19.fa | awk '{print substr($1,2)}' > mapfile

MAPFILE=mapfile

/global/scratch/users/pierrj/eccDNA/pipeline_tests/jithin_issue/git/ecc_caller/generate_bam_file.sh -g ${GENOME_BWA} -1 ${SAMPLE}.sra_1.fastq -2 ${SAMPLE}.sra_2.fastq -s ${SAMPLE} -t ${THREADS} -m ${MAPFILE}


/global/scratch/users/pierrj/eccDNA/pipeline_tests/jithin_issue/git/ecc_caller/call_ecc_regions.sh -m ${MAPFILE} -s ${SAMPLE} -t ${THREADS} \
    -b uniq.filtered.sorted.${SAMPLE}.bam -q multimapped.filtered.name_sorted.${SAMPLE}.bam

/global/scratch/users/pierrj/eccDNA/pipeline_tests/jithin_issue/git/ecc_caller/assign_confidence.sh -m ${MAPFILE} -s ${SAMPLE} -t ${THREADS} \
    -b no_secondary.filtered.sorted.${SAMPLE}.bam -r ${SAMPLE}.confirmedsplitreads.bed