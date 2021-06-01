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

python /global/home/users/pierrj/git/python/pav_quantification_missing_genes_first_pass.py guy11_vs_others ${GENE_GFF} ${GENE_CDS}

while read genome; do
blastn -query ${genome}_missing_genes.fasta -db genomes/${genome}.blastdb \
    -num_threads ${SLURM_NTASKS} -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  stitle " \
    -max_target_seqs 1 | awk '{print $1}' | sort | uniq > not_missing_genes/${genome}.not_missing_genes
done < moryzae_oryza_genomes_mapfile

python /global/home/users/pierrj/git/python/pav_quantification_missing_genes_first_pass.py  guy11_vs_others not_missing_genes ${GENE_GFF} pav_scores_per_gene

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav/G3.all_eccdnas_mapfile

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath ${sample}.ecc_caller_out.splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav

/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -n ${OUTPUT_NAME} \
    -b ${GENE_BEDFILE} -i ${SAMPLE_MAPFILE} -e ${ECCDNA_MAPFILE}