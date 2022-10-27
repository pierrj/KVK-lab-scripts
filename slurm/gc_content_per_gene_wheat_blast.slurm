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

    bedtools nuc -fi ${GENOME_FASTA} -bed ${GENES_BED} | awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' > gc_content_per_gene/${genome}.gene.gc

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