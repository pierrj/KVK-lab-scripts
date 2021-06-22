#!/bin/bash
#SBATCH --job-name=circos_plot
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot

GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed
COPIA_FILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_copia_locs.bed
GYPSY_FILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_gypsy_locs.bed
REPEAT_FILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_all_tes_wrepeats.bed
CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
EFFECTOR_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justeffectors.bed
PAV_SCORES=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav/pav_scores_per_gene
RPKM_FILE=/global/scratch/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/expression_v_eccdnas/micro_dnas/G3.RPKMvsSRs

bedtools makewindows -g ${CHROM_SIZES} -w 1000 > guy11_genome_baoetal2017.1000bpwindows

bedtools makewindows -g ${CHROM_SIZES} -w 10000 > guy11_genome_baoetal2017.10kbpwindows

## tracks to make

# effector locs

input_file=${EFFECTOR_BEDFILE}
track=effector_locs

grep -v '#' ${input_file} | awk -v OFS='\t' '{print substr($1,11,2), $2, $3}' | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > ${track}

# gypsy locs

input_file=${GYPSY_FILE}
track=gypsy_locs

grep -v '#' ${input_file} | awk -v OFS='\t' '{print substr($1,11,2), $2, $3}' | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > ${track}

# copia locs

input_file=${COPIA_FILE}
track=copia_locs

grep -v '#' ${input_file} | awk -v OFS='\t' '{print substr($1,11,2), $2, $3}' | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > ${track}

# PAV per gene
input_file=${PAV_SCORES}
track=pav_per_gene

paste <(awk '{print $1, $2, $3}' ${GENE_BEDFILE}) <(awk '{print $2}' ${input_file}) | awk -v OFS='\t' '{if ($4 == 162) {print substr($1,11,2), $2, $3, 1} else {print substr($1,11,2), $2, $3, 0}}' > ${track}

# expression

input_file=${RPKM_FILE}
track=rpkm

paste <(awk '{print $1, $2, $3}' ${GENE_BEDFILE}) <(awk '{print $2}' ${input_file}) | awk -v OFS='\t' '{print $1, $2, $3, log($4)/log(10)}' > ${track}_old

max=$(awk '{print $4}' ${track}_old | sort -k1,1n | tail -1)

awk -v OFS='\t' '{print $1, $2, $3, $4/max}' ${track}_old > ${track}

# large ecc density

track=large_ecc_density
group=large_eccdnas

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample
do
cd ${sample}
bedtools coverage -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/guy11_genome_baoetal2017.1000bpwindows \
    -b ${group}_splitreads -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.genomecoverage.bed
cd ..
done < mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

while read sample; do
cd ${sample}
total=$(wc -l ${group}_splitreads | awk '{print $1/1000000}')
echo -e ${sample}'\t'${total} >> ../${group}.normalization_table
cd ..
done < mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.${group}.genomecoverage.bed -n ${group}.normalization_table -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

awk -v OFS='\t' '{print $1, $2, log($3)/log(10)}' G3.normalized_binned > /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/${track}_old

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot

max=$(awk '{print $3}' ${track}_old | sort -k1,1n | tail -1)

awk -v OFS='\t' -v max=$max '{OFMT="%f";print substr($1,11,2), $2-1000, $2, $3/max}' ${track}_old | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > ${track}

# micro dna density

track=micro_dna_density
group=micro_dnas

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample
do
cd ${sample}
bedtools coverage -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/guy11_genome_baoetal2017.1000bpwindows \
    -b ${group}_splitreads -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.genomecoverage.bed
cd ..
done < mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

while read sample; do
cd ${sample}
total=$(wc -l ${group}_splitreads | awk '{print $1/1000000}')
echo -e ${sample}'\t'${total} >> ../${group}.normalization_table
cd ..
done < mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.${group}.genomecoverage.bed -n ${group}.normalization_table -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

awk -v OFS='\t' '{print $1, $2, log($3)/log(10)}' G3.normalized_binned > /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/${track}_old

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot

max=$(awk '{print $3}' ${track}_old | sort -k1,1n | tail -1)

awk -v OFS='\t' -v max=$max '{OFMT="%f";print substr($1,11,2), $2-1000, $2, $3/max}' ${track}_old | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > ${track}

# ltr ecc density

track=ltr_ecc_density
group=ltr_eccs

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample
do
cd ${sample}
bedtools coverage -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/guy11_genome_baoetal2017.1000bpwindows \
    -b ${group}_splitreads -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.genomecoverage.bed
cd ..
done < mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

while read sample; do
cd ${sample}
total=$(wc -l ${group}_splitreads | awk '{print $1/1000000}')
echo -e ${sample}'\t'${total} >> ../${group}.normalization_table
cd ..
done < mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.${group}.genomecoverage.bed -n ${group}.normalization_table -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

awk -v OFS='\t' '{print $1, $2, log($3)/log(10)}' G3.normalized_binned > /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/${track}_old

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot

max=$(awk '{print $3}' ${track}_old | sort -k1,1n | tail -1)

awk -v OFS='\t' -v max=$max '{OFMT="%f";print substr($1,11,2), $2-1000, $2, $3/max}' ${track}_old | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > ${track}
