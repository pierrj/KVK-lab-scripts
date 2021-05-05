#!/bin/bash
#SBATCH --job-name=map_pacbio_get_chimeras
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/chimeric_eccs

while read sample; do
    samtools fastq /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all/${sample}/${sample}_pacbio.bam > ${sample}_pacbio.fastq
    /global/scratch/users/pierrj/minimap2/minimap2 -ax asm20 /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_pacbio/guy11_genome_baoetal2017.mmi ${sample}_pacbio.fastq > ${sample}_aligned.bam
    /global/home/users/pierrj/git/bash/get_chimeric_eccs.sh -s ${sample} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/uniq.filtered.sorted.${sample}.bam -p ${sample}_aligned.bam
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all/mapfile