#!/bin/bash
#SBATCH --job-name=eccdnas_guy11
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

genome=guy11

module purge
module load bwa
module load samtools
module load bedtools

PROJ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/${genome}_eccdnas
THREADS=$SLURM_NTASKS

if [ -d "$PROJ_DIR" ]; then
    rm -r "$PROJ_DIR"
fi

mkdir $PROJ_DIR

cd $PROJ_DIR

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
GFF_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/GUY11_fungap_out.fixed.gff3
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

## get gene names and lengths
awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' $GFF_FILE > ${genome}.genes.bed
awk '{print $4}' ${genome}.genes.bed > ${genome}.genes_column
awk '{print $3-$2}' ${genome}.genes.bed > ${genome}.gene_lengths

bwa index ${GENOME_FASTA}

while read sample; do
    # map eccdna reads to reference genome and sort
    bwa mem -t ${THREADS} $GENOME_FASTA /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}_R1.fastq \
    /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}_R2.fastq | \
        samtools sort -@ ${THREADS} -o ${sample}.bam -
    samtools index ${sample}.bam
    ## run bedtools coverage on the mapped histone mark files, normalize by read count
    bedtools coverage -sorted -counts -a ${genome}.genes.bed -g ${GENOME_CHROMSIZES} -b ${sample}.bam | awk '{print $(NF)}' > ${sample}.coverage
    # calculate RPK values for all genes
    paste ${sample}.coverage ${genome}.gene_lengths | awk '{print $1/($2)}' > ${sample}.coverage.RPK
    # sum up RPKs total
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${sample}.coverage.RPK)
    # calculate per million RPKs
    awk -v N=$sum_rpks '{print $1/N}' ${sample}.coverage.RPK > ${sample}.coverage.normalized
done < $MAPFILE

## average all samples together
paste *.coverage.normalized | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > all.coverage
paste ${genome}.genes_column all.coverage > ${genome}_eccdnacov_per_gene.txt