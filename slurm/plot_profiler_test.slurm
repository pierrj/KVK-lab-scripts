#!/bin/bash
#SBATCH --job-name=plot_profiler_test
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/plot_profiler

REGIONS_BED=ERR4832458.dels.bed
DENSITY_BED=GCA_000002495.2_MG8_genomic.fna.out.gff
GENOME=moryzae_70-15_ref_with_mito.fasta
WINDOWS=10
THREADS=${SLURM_NTASKS}
OUTPUT_NAME=test

genome_basename=$(basename ${GENOME})

samtools faidx ${GENOME_FILE}
cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
CHROM_SIZES=${GENOME_FILE}.chromsizes

bedtools makewindows -g ${CHROM_SIZES} -w ${WINDOWS} > ${genome_basename}.${WINDOWS}windows

density_bed_basename=$(basename ${DENSITY_BED} )

bedtools coverage -a ${genome_basename}.${WINDOWS}windows \
    -b ${DENSITY_BED} -g ${CHROM_SIZES} | \
    TEST WHAT THE OUTPUT ACTUALLY LOOKS LIKE
    awk -v OFS='\t' '{print $1, $3, $4}' > ${density_bed_basename}.genomecoverage.bed