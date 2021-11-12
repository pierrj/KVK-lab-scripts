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
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav

# gff file of genes
GENE_GFF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
# fasta file of genes
GENE_CDS=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.genes.fasta
# list of all rice infecting moryzae genomes
GENOMES_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_re_check/moryzae_oryza_genomes_mapfile
# blast db versions of genomes
GENOME_BLASTDB_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_re_check/genomes

# process dir from orthofinder comparing guy11 to other genomes for all orthogroups
# outputs candidate missing genes for each genome and their sequences
python /global/home/users/pierrj/git/python/pav_quantification_missing_genes_first_pass.py guy11_vs_others ${GENE_GFF} ${GENE_CDS}

if [ -d "not_missing_genes" ]; then
    rm -r not_missing_genes
fi

mkdir not_missing_genes

# make sure genes sequences are actually missing from the genome (like you would expect from an eccdna deletion)
# this is very conservative, so any blast hit means the gene isn't actually missing
while read genome; do # loop through genomes
blastn -query ${genome}_missing_genes.fasta -db ${GENOME_BLASTDB_DIR}/${genome}.blastdb \
    -num_threads ${SLURM_NTASKS} -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  stitle " \
    -max_target_seqs 1 | awk '{print $1}' | sort | uniq > not_missing_genes/${genome}.not_missing_genes
done < ${GENOMES_MAPFILE}

# for each orthogroup that are in guy11, count if they are missing in other genomes or not and how many times
python /global/home/users/pierrj/git/python/pav_quantification_missing_genes_second_pass.py  guy11_vs_others not_missing_genes ${GENE_GFF} pav_scores_per_gene

# mapfile containing list of eccdna locations
ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav/G3.all_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

# just copy all junction split read files
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

# output never (eccdna-absent) and common (eccdna-associated) genes using encompassing split reads
/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -n ${OUTPUT_NAME} \
    -b ${GENE_BEDFILE} -i ${SAMPLE_MAPFILE} -e ${ECCDNA_MAPFILE}