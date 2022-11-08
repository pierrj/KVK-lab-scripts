#!/bin/bash
#SBATCH --job-name=gc_gene_stats
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

DATA_OUTPUT_DIR=gc_content

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

source activate /global/scratch/users/pierrj/conda_envs/deeptools

bedtools nuc -fi ${GENOME_FASTA} -bed ${PROTEIN_BED} | awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' > gene.gc

bedtools flank -g ${GENOME_CHROM_SIZES} -i ${PROTEIN_BED} -b 1000 | \
    bedtools nuc -fi ${GENOME_FASTA} -bed - | \
    awk -v OFS='\t' '{ if (NR > 1) {print $4, $5}}' | \
    awk -v OFS='\t' '{ sum[$1] += $2; count[$1] += 1 } END { for ( key in count ) { print key, sum[key] / count[key] } }' > 1kb.gc

cut -f1 1kb.gc > gene_col

cut -f2 gene.gc > gene.gc_col
cut -f2 1kb.gc > 1kb.gc_col

paste gene_col gene.gc_col 1kb.gc_col > gc_table.txt