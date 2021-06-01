#!/bin/bash
#SBATCH --job-name=eccdna_pav_scores
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav

GENE_GFF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
GENE_CDS=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.genes.fasta
GENOMES_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_re_check/moryzae_oryza_genomes_mapfile
GENOME_BLASTDB_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_re_check/genomes

python /global/home/users/pierrj/git/python/pav_quantification_missing_genes_first_pass.py guy11_vs_others ${GENE_GFF} ${GENE_CDS}

if [ -d "not_missing_genes" ]; then
    rm -r not_missing_genes
fi

mkdir not_missing_genes

while read genome; do
blastn -query ${genome}_missing_genes.fasta -db ${GENOME_BLASTDB_DIR}/${genome}.blastdb \
    -num_threads ${SLURM_NTASKS} -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  stitle " \
    -max_target_seqs 1 | awk '{print $1}' | sort | uniq > not_missing_genes/${genome}.not_missing_genes
done < ${GENOMES_MAPFILE}

python /global/home/users/pierrj/git/python/pav_quantification_missing_genes_second_pass.py  guy11_vs_others not_missing_genes ${GENE_GFF} pav_scores_per_gene

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav/G3.all_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath ${sample}.ecc_caller_out.splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav

OUTPUT_NAME=G3
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile

/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -n ${OUTPUT_NAME} \
    -b ${GENE_BEDFILE} -i ${SAMPLE_MAPFILE} -e ${ECCDNA_MAPFILE}