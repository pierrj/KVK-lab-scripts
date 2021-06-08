#!/bin/bash
#SBATCH --job-name=eccdnavsrnaseqcomparison_ct_if
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
SAMPLE=IF
SRA_LIST=/global/scratch/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/zhang_et_al_rnaseq_accessions.txt
THREADS=${SLURM_NTASKS}
LIB_TYPE=1
GFF_FILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
SAMPLE_MAPFILE=/global/scratch/pierrj/eccDNA/magnaporthe_pureculture/illumina_if_v_g3/sample_mapfile_IF
copia_file=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_copia_locs.bed
gypsy_file=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_gypsy_locs.bed

module load STAR





ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas_ct_if/large_eccdnas/G3.large_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina_if_v_g3

while read sample; do
cd ${sample}
    realpath large_eccdnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile_if

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas_ct_if/large_eccdnas/

/global/home/users/pierrj/git/bash/ecc_dna_vs_expression.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n a





ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas_ct_if/micro_dnas/G3.micro_dnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina_if_v_g3

while read sample; do
cd ${sample}
    realpath micro_dnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile_if

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas_ct_if/micro_dnas/

/global/home/users/pierrj/git/bash/ecc_dna_vs_expression.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n a







ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas_ct_if/gene_containing_eccdnas/G3.all_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina_if_v_g3

while read sample; do
cd ${sample}
    realpath ${sample}.ecc_caller_out.splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile_if

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas_ct_if/gene_containing_eccdnas

/global/home/users/pierrj/git/bash/ecc_dna_vs_expression.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n n