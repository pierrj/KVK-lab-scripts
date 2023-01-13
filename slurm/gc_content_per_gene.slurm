#!/bin/bash
#SBATCH --job-name=gc_content_per_gene
#SBATCH --partition=savio2
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

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/

source activate /global/scratch/users/pierrj/conda_envs/deeptools

cp -r /global/scratch/users/pierrj/fungap_runs/gladieux_all/all_gffs_fixed .

if [ ! -d all_beds_fixed ]; then
    mkdir all_beds_fixed
fi

for file in all_gffs_fixed/*; do
    genome=$(basename $file | awk '{print substr($1, 0, length($1)-22)}')
    echo $genome
    awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' $file > all_beds_fixed/${genome}.genes.bed
done

# while read genome; do
#     samtools faidx ${genome}.fasta
#     cut -f1-2 ${genome}.fasta.fai > ${genome}.chromsizes
# done < ../assemblies_list


if [ ! -d gc_content_per_gene ]; then
    mkdir gc_content_per_gene
fi

if [ ! -d gc_content_1kb_flank ]; then
    mkdir gc_content_1kb_flank
fi

if [ ! -d gc_content_table ]; then
    mkdir gc_content_table
fi

while read genome; do
    echo ${genome}
    GENOME_FASTA=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/assemblies/${genome}.fasta
    GENOME_CHROM_SIZES=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/assemblies/${genome}.chromsizes
    GENES_BED=all_beds_fixed/${genome}.genes.bed

    # get gc content per gene
    bedtools nuc -fi ${GENOME_FASTA} -bed ${GENES_BED} > gc_content_per_gene/${genome}.gene.gc

    # get 1kb flanking per gene, on either side
    bedtools flank -g ${GENOME_CHROM_SIZES} -i ${GENES_BED} -b 1000 | \
        bedtools nuc -fi ${GENOME_FASTA} -bed - > gc_content_1kb_flank/${genome}.1kb.gc
    
    ## need to manually calculate gc content here to deal with Ns
    ## also need to deal with flanking regions that are too small or have too many Ns
    ## basically filtering regions less than 100 bp and regions that have less than 10% Ns
    python /global/home/users/pierrj/git/python/generate_gc_table_per_genome.py gc_content_per_gene/${genome}.gene.gc gc_content_1kb_flank/${genome}.1kb.gc $genome
    
done < genomes_mapfile_no_mgrisae

cat gc_content_table/* > all_genomes.gc_table.txt