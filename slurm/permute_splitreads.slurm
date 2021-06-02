#!/bin/bash
#SBATCH --job-name=permute_splitreads
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


LTR_TE_LOCATIONS=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile
OUTPUT_NAME=G3

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/permuted_splitreads/G3.nonltr_splitreads.mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    cat micro_dnas_splitreads large_eccdnas_splitreads > /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/permuted_splitreads/${sample}.nonltr_splitreads.bed
    echo ${sample}.nonltr_splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/permuted_splitreads

for i in $(seq 0 999); do
    echo ${i}
    SHUFFLED_ECCDNA_MAPFILE=shuffled_eccdnas_${i}
    if [ -f "${SHUFFLED_ECCDNA_MAPFILE}" ]; then
        rm ${SHUFFLED_ECCDNA_MAPFILE}
    fi
    while read ecc_file; do
        bedtools shuffle -i ${ecc_file} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > ${ecc_file}_shuffled_${i}
        echo ${ecc_file}_shuffled_${i} >> ${SHUFFLED_ECCDNA_MAPFILE}
    done < ${ECCDNA_MAPFILE}
    /global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -n ${OUTPUT_NAME} \
        -b ${GENE_BEDFILE} -i ${SAMPLE_MAPFILE} -e ${SHUFFLED_ECCDNA_MAPFILE}
    awk '{print $2}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${i}.shuffled.splitreadspergene
    cp ${OUTPUT_NAME}.common.genes ${i}.shuffle.common.genes
done

paste ${OUTPUT_NAME}.allgenenames $(find . -maxdepth 1 -name "*.shuffled.splitreadspergene" | xargs -r ls -1 | tr "\n" " ") > shuffled.splitreadspergene.fulltable

paste $(find . -maxdepth 1 -name "*.shuffle.common.genes" | xargs -r ls -1 | tr "\n" " ") > shuffled.common.genes.fulltable