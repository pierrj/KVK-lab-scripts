#!/bin/bash
#SBATCH --job-name=gc_content_per_gene_wheatblast
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

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest

MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/wheat_blast_busco_greater_than_90_annotated

source activate /global/scratch/users/pierrj/conda_envs/deeptools

cp -r /global/scratch/users/pierrj/fungap_runs/wheat_blast/all_gffs_fixed .

if [ ! -d all_beds_fixed ]; then
    mkdir all_beds_fixed
fi

for file in all_gffs_fixed/*; do
    genome=$(basename $file | awk '{print substr($1, 0, length($1)-22)}')
    echo $genome
    awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' $file > all_beds_fixed/${genome}.genes.bed
done

while read genome; do
    echo ${genome}
    GENOME_FASTA=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/${genome}_genomic.fna
    samtools faidx ${GENOME_FASTA}
    cut -f1-2 ${GENOME_FASTA}.fai > ${GENOME_FASTA}.chromsizes
done < ${MAPFILE}


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
    GENOME_FASTA=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/${genome}_genomic.fna
    GENOME_CHROM_SIZES=${GENOME_FASTA}.chromsizes
    GENES_BED=all_beds_fixed/${genome}.genes.bed

    # get gc content per gene
    ## THESE ACTUALLY CALCULATE AT CONTENT NOT GC CONTENT ##
    bedtools nuc -fi ${GENOME_FASTA} -bed ${GENES_BED} | awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' > gc_content_per_gene/${genome}.gene.gc


    # get 1kb flanking per gene, on either side
    # then get region gc content
    # then average the two regions with awk
    bedtools flank -g ${GENOME_CHROM_SIZES} -i ${GENES_BED} -b 1000 | \
        bedtools nuc -fi ${GENOME_FASTA} -bed - | \
        awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' | \
        awk -v OFS='\t' '{ sum[$1] += $2; count[$1] += 1 } END { for ( key in count ) { print key, sum[key] / count[key] } }' > gc_content_1kb_flank/${genome}.1kb.gc

    cut -f1 gc_content_1kb_flank/${genome}.1kb.gc > ${genome}.gene_col

    cut -f2 gc_content_per_gene/${genome}.gene.gc > ${genome}.gene.gc_col
    cut -f2 gc_content_1kb_flank/${genome}.1kb.gc > ${genome}.1kb.gc_col

    paste ${genome}.gene_col ${genome}.gene.gc_col ${genome}.1kb.gc_col > gc_content_table/${genome}.gc_table.txt
    
    rm ${genome}*
done < ${MAPFILE}

cat gc_content_table/* > all_genomes.gc_table.txt