#!/bin/bash
#SBATCH --job-name=eccDNA_gene_analysis_wrapper
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
source activate pfamscan

module remove r

module load r/3.4.2


cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_postmapq0/G3_testrun

rm -r *

cp /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile .
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile .

while read sample; do
cp /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.confirmedsplitreads.bed .
done < mapfile

PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_output_fungap_guy11_12_28_20.txt
OUTPUT_NAME=G3
SCORE_CUTOFF=0.6
MAPFILE=mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
TE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_pacbio_ET_KROJ.features.justtranpsons.bed
PFAM_DIR=/global/scratch/users/pierrj/pfam/
CDS_FASTA=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20_prot.faa
THREADS=${SLURM_NTASKS}
SAMPLE_MAPFILE=G3.sample_mapfile

/global/home/users/pierrj/git/bash/gene_analysis_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME} -c ${SCORE_CUTOFF} -m ${MAPFILE} -b ${GENE_BEDFILE} -e ${TE_BEDFILE} \
            -f ${PFAM_DIR} -a ${CDS_FASTA} -t ${THREADS} -i ${SAMPLE_MAPFILE}