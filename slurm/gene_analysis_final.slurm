#!/bin/bash
#SBATCH --job-name=gene_analysis_final
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis

source activate pfamscan

module remove r

module load r/3.4.2

PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_output_fungap_guy11_12_28_20.txt
OUTPUT_NAME=G3
SCORE_CUTOFF=0.6
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
TE_BEDFILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_all_tes_norepeats.bed
REPEAT_BEDFILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_all_tes_wrepeats.bed
PFAM_DIR=/global/scratch/users/pierrj/pfam/
CDS_FASTA=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20_prot.faa
THREADS=${SLURM_NTASKS}
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile


ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/G3.all_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath ${sample}.ecc_caller_out.splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis

/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -n ${OUTPUT_NAME} \
    -b ${GENE_BEDFILE} -i ${SAMPLE_MAPFILE} -e ${ECCDNA_MAPFILE}


cp ${OUTPUT_NAME}.common.genes common/${OUTPUT_NAME}.common.genes

cp ${OUTPUT_NAME}.neverfound.genes neverfound/${OUTPUT_NAME}.neverfound.genes

cp ${OUTPUT_NAME}.allgenenames common/${OUTPUT_NAME}.allgenenames

cp ${OUTPUT_NAME}.allgenenames neverfound/${OUTPUT_NAME}.allgenenames

cp ${GENE_BEDFILE} common/
cp ${GENE_BEDFILE} neverfound/

cp ${TE_BEDFILE} common/
cp ${TE_BEDFILE} neverfound/

cp ${REPEAT_BEDFILE} common/
cp ${REPEAT_BEDFILE} neverfound/

cp ${PANNZER_INPUT} common/
cp ${PANNZER_INPUT} neverfound/


python /global/home/users/pierrj/git/python/parse_pannzer_output.py ${PANNZER_INPUT} ${SCORE_CUTOFF} ${OUTPUT_NAME}.filtered.GO.txt

awk -v OFS='\t' '{print $0}' ${OUTPUT_NAME}.filtered.GO.txt > ${OUTPUT_NAME}.filtered.GO.awk.txt

cp ${OUTPUT_NAME}.filtered.GO.awk.txt common/
cp ${OUTPUT_NAME}.filtered.GO.awk.txt neverfound/


cd common

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${OUTPUT_NAME}.common.genes -a ${OUTPUT_NAME}.allgenenames -n ${OUTPUT_NAME}.common_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

cd ..


cd neverfound

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${OUTPUT_NAME}.neverfound.genes -a ${OUTPUT_NAME}.allgenenames -n ${OUTPUT_NAME}.never_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

cd ..