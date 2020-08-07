#!/bin/bash
#SBATCH --job-name=human_map_time_check
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd ${SLURM_SUBMIT_DIR}

mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample="SRR6315407"
bamfile="/global/scratch/users/pierrj/eccDNA/2018_moller/SRR6315407/SRR6315407.sorted.mergedandpe.bwamem.bam"

cd ${sample}

date +"%T"

/global/home/users/pierrj/git/bash/generate_bam_file.sh -g ${genome_bwa} \
 -1 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_R1.fastq -2 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_R2.fastq \
 -s ${sample} -t ${SLURM_NTASKS} -m ${mapfile}

date +"%T"

/global/home/users/pierrj/git/bash/call_ecc_regions.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b filtered.sorted.${sample}.bam

date +"%T"

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b filtered.sorted.${sample}.bam

date +"%T"