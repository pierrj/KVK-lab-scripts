#!/bin/bash
#SBATCH --job-name=is_pyret_active
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/is_pyret_active/

genome_bwa='/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa'
sample=G3_1A
mapfile='/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames'
bed='/global/scratch/users/pierrj/repeatmasker/moryzae_ltr_parts/pyret_locs.bed'

# /global/home/users/pierrj/git/bash/generate_bam_file.sh -g ${genome_bwa} \
#  -1 /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R1.fastq -2 /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A_R2.fastq \
#  -s ${sample} -t ${SLURM_NTASKS} -m ${mapfile}

# samtools depth -a filtered.sorted.${sample}.bam > ${sample}.genomecoverage

# /global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t G3_1A.genomecoverage -n filtered.sorted.G3_1A.bam

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1000000 -b 100

cd G3_1A

/global/home/users/pierrj/git/bash/generate_coverageplot_withhighlighting.sh -m ${mapfile} -s ${sample} -c ${sample}.normalized_binned  -b ${bed}

date +"%T"

bedtools genomecov -d -ibam filtered.sorted.${sample}.bam > ${SAMPLE}.genomecoverage.bedtools

date +"%T"

samtools depth -a filtered.sorted.${sample}.bam > ${sample}.genomecoverage

date +"%T"