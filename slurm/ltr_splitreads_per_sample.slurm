#!/bin/bash
#SBATCH --job-name=ltr_splitreads_per_sample
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads

TE_BEDFILE=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_nolengthfilter_individual_elements.bed
ELEMENT_MAPFILE=/global/scratch/users/pierrj/repeatmasker/moryzae_final/mapfile_elements

while read sample
do
bamfile=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam
/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${sample} -b ${bamfile}
split_read_file=candidate_ltr_srs.${sample}.bam
/global/home/users/pierrj/git/bash/get_ltr_sr_circle_count_all_elements.sh -b ${TE_BEDFILE} \
    -m ${ELEMENT_MAPFILE} \
    -s ${split_read_file} -o ${sample}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

## normalize and average

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

while read sample; do
    if [ -d "${sample}" ]; then
        rm -r ${sample}
    fi
    mkdir ${sample}
    mv ${sample}.ltr_sr_cov_perfeature ${sample}/${sample}.ltr_sr_cov_perfeature
    mv ${sample}.read_cov_perfeature ${sample}/${sample}.read_cov_perfeature
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

if [ -f "G3.normalization_table.ltr_sr_cov_perfeature" ]; then
    rm G3.normalization_table.ltr_sr_cov_perfeature
fi

while read sample; do
    total=$( awk '{sum += $2} END {print sum}' ${sample}/${sample}.ltr_sr_cov_perfeature | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${total} >> G3.normalization_table.ltr_sr_cov_perfeature
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.ltr_sr_cov_perfeature -n G3.normalization_table.ltr_sr_cov_perfeature -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n

mv G3.normalized_binned ltr_sr_cov_perfeature_normalized


if [ -f "G3.normalization_table.read_cov_perfeature" ]; then
    rm G3.normalization_table.read_cov_perfeature
fi

while read sample; do
    total=$( awk '{sum += $2} END {print sum}' ${sample}/${sample}.read_cov_perfeature | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${total} >> G3.normalization_table.read_cov_perfeature
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.ltr_sr_cov_perfeature -n G3.normalization_table.read_cov_perfeature -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n

mv G3.normalized_binned read_cov_perfeature_normalized


paste ltr_sr_cov_perfeature_normalized read_cov_perfeature_normalized > ltr_sr_and_read_cov_perfeature_normalized