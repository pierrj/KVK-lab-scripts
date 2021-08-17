#!/bin/bash
#SBATCH --job-name=uniq_eccs_moryzae
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads/uniq_eccs

SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"

while read sample; do
    if [ -d "${sample}" ]; then
        rm -r ${sample}
    fi
    mkdir ${sample}
    cd ${sample}
    /global/home/users/pierrj/git/bash/ecc_caller_only_unique.sh -m ${mapfile} \
        -s ${sample} \
        -t ${SLURM_NTASKS} \
        -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/uniq.filtered.sorted.${sample}.bam
    /global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
        -s ${sample} \
        -t ${SLURM_NTASKS} \
        -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam \
        -r ${sample}.unique.confirmedsplitreads.bed
    cp ${sample}.ecc_caller_out.details.txt ${sample}.ecc_caller_out.uniq.details.txt 
    cd ..
done < ${SAMPLE_MAPFILE}