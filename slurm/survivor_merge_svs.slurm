#!/bin/bash
#SBATCH --job-name=survivor_merge_svs
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/survivor

REFERENCE=moryzae_70-15_ref_with_mito.fasta
SAMPLE=ERR4757126
GENOME_DB=moryzae_70-15_ref_with_mito.fasta
READONE=${SAMPLE}.sra_1.fastq
READTWO=${SAMPLE}.sra_2.fastq

if [ -f ${SAMPLE}.survivor.vcflist ]; then
    rm ${SAMPLE}.survivor.vcflist
fi

cp 

realpath /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/delly/${SAMPLE}.delly.vcf >> ${SAMPLE}.survivor.vcflist
realpath /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/lumpy/${SAMPLE}.lumpy.vcf >> ${SAMPLE}.survivor.vcflist
realpath /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/manta/${SAMPLE}.manta.vcf >> ${SAMPLE}.survivor.vcflist
realpath /global/scratch/users/pierrj/sv_calling_moryzae/my_tests/wham/${SAMPLE}.wham.vcf >> ${SAMPLE}.survivor.vcflist


## in order
# distance between breakpoints
# number of supporting callers
# take type into account
# take strand into account
# disabled?
# minimum size
# output

/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge ${SAMPLE}.survivor.vcflist 1000 3 1 0 0 50 ${SAMPLE}.all.vcf

some more filtering needed here such as SV length and mapping quality of regions