#!/bin/bash
#SBATCH --job-name=eccs_per_ltr_per_sample
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/

LTR_TE_FILE=/global/scratch/users/pierrj/repeatmasker/moryzae_final/copia_int_renamed.only_LTRs_of_interest.identitycutoff.bed
ELEMENT_MAPFILE=/global/scratch/users/pierrj/repeatmasker/moryzae_final/mapfile_elements

while read SAMPLE
do
cd ${SAMPLE}
    if [ -f "${SAMPLE}.sr_count_per_element" ]; then
        rm ${SAMPLE}.sr_count_per_element
    fi
    while read ELEMENT
    do
        grep ${ELEMENT} ${LTR_TE_FILE} > ${ELEMENT}.loc.bed
        bedtools intersect -wao -a ${SAMPLE}.confirmedsplitreads.bed -b ${ELEMENT}.loc.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > ${SAMPLE}.${ELEMENT}.ltr_eccs
        sr_count=$(wc -l ${SAMPLE}.${ELEMENT}.ltr_eccs | awk '{print $1}')
        echo -e ${ELEMENT}'\t'${sr_count} >> ${SAMPLE}.sr_count_per_element
    done < ${ELEMENT_MAPFILE}
cd ..
done < mapfile_IF

if [ -f "IF.normalization_table.ltr_splitreads" ]; then
    rm IF.normalization_table.ltr_splitreads
fi

while read sample; do
cd ${sample}
sum=$(wc -l ${sample}.confirmedsplitreads.bed | awk '{print $1}')
echo -e ${sample}'\t'${sum} >> ../IF.normalization_table.ltr_splitreads
cd ..
done < mapfile_IF

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m mapfile_IF

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t IF_1A.sr_count_per_element -n IF.normalization_table.ltr_splitreads -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1000000 -b 1 -c 2 -n n




cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read SAMPLE
do
cd ${SAMPLE}
    if [ -f "${SAMPLE}.sr_count_per_element" ]; then
        rm ${SAMPLE}.sr_count_per_element
    fi
    while read ELEMENT
    do
        grep ${ELEMENT} ${LTR_TE_FILE} > ${ELEMENT}.loc.bed
        bedtools intersect -wao -a ${SAMPLE}.confirmedsplitreads.bed -b ${ELEMENT}.loc.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > ${SAMPLE}.${ELEMENT}.ltr_eccs
        sr_count=$(wc -l ${SAMPLE}.${ELEMENT}.ltr_eccs | awk '{print $1}')
        echo -e ${ELEMENT}'\t'${sr_count} >> ${SAMPLE}.sr_count_per_element
    done < ${ELEMENT_MAPFILE}
cd ..
done < mapfile

if [ -f "G3.normalization_table.ltr_splitreads" ]; then
    rm G3.normalization_table.ltr_splitreads
fi

while read sample; do
cd ${sample}
sum=$(wc -l ${sample}.confirmedsplitreads.bed | awk '{print $1}')
echo -e ${sample}'\t'${sum} >> ../G3.normalization_table.ltr_splitreads
cd ..
done < mapfile

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t G3_1A.sr_count_per_element -n G3.normalization_table.ltr_splitreads -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1000000 -b 1 -c 2 -n n