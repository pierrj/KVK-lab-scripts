#!/bin/bash
#SBATCH --job-name=rarefaction_G3_1A
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

SAMPLE=G3_1A


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/rarefaction_tech_reps/${SAMPLE}

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"

readcount=$(wc -l /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/${SAMPLE}_R1.fastq)

for i in $(seq 0.1 0.1 1.0); do

    sample=${i}

    if [ -d "${i}" ]; then
        rm -r ${i}
    fi

    mkdir ${i}

    cd ${i}

    samtools view -s ${i} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/uniq.filtered.sorted.${SAMPLE}.bam > uniq.filtered.sorted.${sample}.bam
    samtools view -s ${i} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/multimapped.filtered.name_sorted.${SAMPLE}.bam > multimapped.filtered.name_sorted.${sample}.bam
    samtools view -s ${i} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/no_secondary.filtered.sorted.${SAMPLE}.bam > no_secondary.filtered.sorted.${sample}.bam

    ## call putative eccDNA forming regions, using both uniquely mapped and multi-mapped reads
    /global/home/users/pierrj/git/bash/ecc_caller_mapq0.sh -m ${mapfile} \
        -s ${sample} \
        -t ${SLURM_NTASKS} \
        -b uniq.filtered.sorted.${sample}.bam \
        -q multimapped.filtered.name_sorted.${sample}.bam

    ## assign confidence for eccDNA regions
    /global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
        -s ${sample} \
        -t ${SLURM_NTASKS} \
        -b no_secondary.filtered.sorted.${sample}.bam \
        -r ${sample}.confirmedsplitreads.bed

    cd ..
done 