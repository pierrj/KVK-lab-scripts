#!/bin/bash
#SBATCH --job-name=gc_content_per_gene_wheat_blast
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/

source activate /global/scratch/users/pierrj/conda_envs/deeptools

MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_annotated_mapfile
GENOMES_PATH=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate

if [ ! -d all_beds ]; then
    mkdir all_beds
fi

if [ ! -d gc_content_per_gene ]; then
    mkdir gc_content_per_gene
fi

if [ ! -d gc_content_1kb_flank ]; then
    mkdir gc_content_1kb_flank
fi

if [ ! -d gc_content_5kb_flank ]; then
    mkdir gc_content_5kb_flank
fi

if [ ! -d gc_content_table ]; then
    mkdir gc_content_table
fi


while read genome; do
    GENOME_FASTA=${GENOMES_PATH}/${genome}_genomic.fna
    GENOME_CHROM_SIZES=${GENOMES_PATH}/${genome}.chromsizes
    GENES_BED=all_beds/${genome}.genes.bed

    awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' /global/scratch/users/pierrj/fungap_runs/wheat_blast/all_gffs_fixed/${genome}_fungap_out.fixed.gff3 > $GENES_BED

    samtools faidx $GENOME_FASTA
    cut -f1-2 $GENOME_FASTA.fai > $GENOME_CHROM_SIZES

    bedtools nuc -fi ${GENOME_FASTA} -bed ${GENES_BED} | awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' > gc_content_per_gene/${genome}.gene.gc

    bedtools flank -g ${GENOME_CHROM_SIZES} -i ${GENES_BED} -b 1000 | \
        bedtools nuc -fi ${GENOME_FASTA} -bed - | \
        awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' | \
        awk -v OFS='\t' '{ sum[$1] += $2; count[$1] += 1 } END { for ( key in count ) { print key, sum[key] / count[key] } }' > gc_content_1kb_flank/${genome}.1kb.gc

    bedtools flank -g ${GENOME_CHROM_SIZES} -i ${GENES_BED} -b 5000 | \
        bedtools nuc -fi ${GENOME_FASTA} -bed - | \
        awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' | \
        awk -v OFS='\t' '{ sum[$1] += $2; count[$1] += 1 } END { for ( key in count ) { print key, sum[key] / count[key] } }' > gc_content_5kb_flank/${genome}.5kb.gc

    cut -f1 gc_content_5kb_flank/${genome}.5kb.gc > ${genome}.gene_col

    cut -f2 gc_content_per_gene/${genome}.gene.gc > ${genome}.gene.gc_col
    cut -f2 gc_content_1kb_flank/${genome}.1kb.gc > ${genome}.1kb.gc_col
    cut -f2 gc_content_5kb_flank/${genome}.5kb.gc > ${genome}.5kb.gc_col

    paste ${genome}.gene_col ${genome}.gene.gc_col ${genome}.1kb.gc_col ${genome}.5kb.gc_col > gc_content_table/${genome}.gc_table.txt
    
    rm ${genome}*
    
done < $MAPFILE

cat gc_content_table/* > all_genomes.gc_table.txt