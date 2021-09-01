#!/bin/bash
#SBATCH --job-name=blast_eccDNA_ends_yeast
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/homology_at_eccdna_ends/yeast


echo 'this doesnt include the mitochondria filtering scheme from before you clown, its only there bc you put it in manually'

if [ -f "large_eccdnas_all" ]; then
    rm large_eccdnas_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/${sample}/large_eccdnas >> large_eccdnas_all; done < /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/mapfile
sort -k1,1 -k2,2n large_eccdnas_all | uniq > large_eccdnas_all_uniq

if [ -f "micro_dnas_all" ]; then
    rm micro_dnas_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/${sample}/micro_dnas >> micro_dnas_all; done < /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/mapfile
sort -k1,1 -k2,2n micro_dnas_all | uniq > micro_dnas_all_uniq

cat large_eccdnas_all_uniq micro_dnas_all_uniq > all_uniq


# INPUT_FILE=large_eccdnas_all_uniq
# OUTPUT_TABLE=large_eccdnas_observed_homology
# GENOME=/global/scratch/users/pierrj/references/Scer_S288C.fasta
# GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/Scer_S288C.fasta.filtered.chromsizes
# LTR_TE_LOCATIONS=/global/scratch/users/pierrj/references/te_annotations/yeast/yeast_copia_and_gypsy_locs.bed.filtered


# if [ -f "${OUTPUT_TABLE}" ]; then
#     rm ${OUTPUT_TABLE}
# fi

# ## deal with being too close to end
# awk -v OFS='\t' '{print $1, $2-100, $2}' ${GENOME_CHROMSIZES} > guy11_genome_baoetal2017.chromends
# bedtools intersect -v -wa -a ${INPUT_FILE} -b guy11_genome_baoetal2017.chromends > input_file_not_close_to_end

# ## deal with being too close to start
# awk -v OFS='\t' '{ if ($2-100>0) {print $1, $2-100, $2+100}}' input_file_not_close_to_end > eccdna_starts
# awk -v OFS='\t' '{ if ($2-100>0) {print $1, $3-100, $3+100}}' input_file_not_close_to_end > eccdna_ends

# bedtools getfasta -fi ${GENOME} -fo eccdna_starts.fasta -bed eccdna_starts
# bedtools getfasta -fi ${GENOME} -fo eccdna_ends.fasta -bed eccdna_ends

# awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
# awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
# paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq

# parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}

# for i in $(seq 0 99); do
#     echo ${i}
#     bedtools shuffle -i ${INPUT_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > ${INPUT_FILE}_shuffle_${i}
#     OUTPUT_TABLE=${INPUT_FILE}_shuffle_${i}_outputtable
#     awk -v OFS='\t' '{print $1, $2-100, $2}' ${GENOME_CHROMSIZES} > guy11_genome_baoetal2017.chromends
#     bedtools intersect -v -wa -a ${INPUT_FILE}_shuffle_${i} -b guy11_genome_baoetal2017.chromends > input_file_not_close_to_end
#     awk -v OFS='\t' '{ if ($2-100>0) {print $1, $2-100, $2+100}}' input_file_not_close_to_end > eccdna_starts
#     awk -v OFS='\t' '{ if ($2-100>0) {print $1, $3-100, $3+100}}' input_file_not_close_to_end > eccdna_ends
#     bedtools getfasta -fi ${GENOME} -fo eccdna_starts.fasta -bed eccdna_starts
#     bedtools getfasta -fi ${GENOME} -fo eccdna_ends.fasta -bed eccdna_ends
#     awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
#     awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
#     paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq
#     parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}
# done

# paste $(find . -maxdepth 1 -name "${INPUT_FILE}_shuffle_*_outputtable" | xargs -r ls -1 | tr "\n" " ") > ${INPUT_FILE}_shuffle_outputtable


# INPUT_FILE=micro_dnas_all_uniq
# OUTPUT_TABLE=micro_dnas_observed_homology
# GENOME=/global/scratch/users/pierrj/references/Scer_S288C.fasta
# GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/Scer_S288C.fasta.filtered.chromsizes
# LTR_TE_LOCATIONS=/global/scratch/users/pierrj/references/te_annotations/yeast/yeast_copia_and_gypsy_locs.bed.filtered

# if [ -f "${OUTPUT_TABLE}" ]; then
#     rm ${OUTPUT_TABLE}
# fi

# ## deal with being too close to end
# awk -v OFS='\t' '{print $1, $2-100, $2}' ${GENOME_CHROMSIZES} > guy11_genome_baoetal2017.chromends
# bedtools intersect -v -wa -a ${INPUT_FILE} -b guy11_genome_baoetal2017.chromends > input_file_not_close_to_end

# ## deal with being too close to start
# awk -v OFS='\t' '{ if ($2-100>0) {print $1, $2-100, $2+100}}' input_file_not_close_to_end > eccdna_starts
# awk -v OFS='\t' '{ if ($2-100>0) {print $1, $3-100, $3+100}}' input_file_not_close_to_end > eccdna_ends

# bedtools getfasta -fi ${GENOME} -fo eccdna_starts.fasta -bed eccdna_starts
# bedtools getfasta -fi ${GENOME} -fo eccdna_ends.fasta -bed eccdna_ends

# awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
# awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
# paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq

# parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}

# for i in $(seq 0 99); do
#     echo ${i}
#     bedtools shuffle -i ${INPUT_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > ${INPUT_FILE}_shuffle_${i}
#     OUTPUT_TABLE=${INPUT_FILE}_shuffle_${i}_outputtable
#     awk -v OFS='\t' '{print $1, $2-100, $2}' ${GENOME_CHROMSIZES} > guy11_genome_baoetal2017.chromends
#     bedtools intersect -v -wa -a ${INPUT_FILE}_shuffle_${i} -b guy11_genome_baoetal2017.chromends > input_file_not_close_to_end
#     awk -v OFS='\t' '{ if ($2-100>0) {print $1, $2-100, $2+100}}' input_file_not_close_to_end > eccdna_starts
#     awk -v OFS='\t' '{ if ($2-100>0) {print $1, $3-100, $3+100}}' input_file_not_close_to_end > eccdna_ends
#     bedtools getfasta -fi ${GENOME} -fo eccdna_starts.fasta -bed eccdna_starts
#     bedtools getfasta -fi ${GENOME} -fo eccdna_ends.fasta -bed eccdna_ends
#     awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
#     awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
#     paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq
#     parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}
# done

# paste $(find . -maxdepth 1 -name "${INPUT_FILE}_shuffle_*_outputtable" | xargs -r ls -1 | tr "\n" " ") > ${INPUT_FILE}_shuffle_outputtable


INPUT_FILE=all_uniq
OUTPUT_TABLE=all_observed_homology
GENOME=/global/scratch/users/pierrj/references/Scer_S288C.fasta
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/Scer_S288C.fasta.filtered.chromsizes
LTR_TE_LOCATIONS=/global/scratch/users/pierrj/references/te_annotations/yeast/yeast_copia_and_gypsy_locs.bed.filtered

if [ -f "${OUTPUT_TABLE}" ]; then
    rm ${OUTPUT_TABLE}
fi

## deal with being too close to end
awk -v OFS='\t' '{print $1, $2-10, $2}' ${GENOME_CHROMSIZES} > guy11_genome_baoetal2017.chromends
bedtools intersect -v -wa -a ${INPUT_FILE} -b guy11_genome_baoetal2017.chromends > input_file_not_close_to_end

## deal with being too close to start
awk -v OFS='\t' '{ if ($2-10>0) {print $1, $2-10, $2+10}}' input_file_not_close_to_end > eccdna_starts
awk -v OFS='\t' '{ if ($2-10>0) {print $1, $3-10, $3+10}}' input_file_not_close_to_end > eccdna_ends

bedtools getfasta -fi ${GENOME} -fo eccdna_starts.fasta -bed eccdna_starts
bedtools getfasta -fi ${GENOME} -fo eccdna_ends.fasta -bed eccdna_ends

awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq

parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}

for i in $(seq 0 99); do
    echo ${i}
    bedtools shuffle -i ${INPUT_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > ${INPUT_FILE}_shuffle_${i}
    OUTPUT_TABLE=${INPUT_FILE}_shuffle_${i}_outputtable
    awk -v OFS='\t' '{print $1, $2-10, $2}' ${GENOME_CHROMSIZES} > guy11_genome_baoetal2017.chromends
    bedtools intersect -v -wa -a ${INPUT_FILE}_shuffle_${i} -b guy11_genome_baoetal2017.chromends > input_file_not_close_to_end
    awk -v OFS='\t' '{ if ($2-10>0) {print $1, $2-10, $2+10}}' input_file_not_close_to_end > eccdna_starts
    awk -v OFS='\t' '{ if ($2-10>0) {print $1, $3-10, $3+10}}' input_file_not_close_to_end > eccdna_ends
    bedtools getfasta -fi ${GENOME} -fo eccdna_starts.fasta -bed eccdna_starts
    bedtools getfasta -fi ${GENOME} -fo eccdna_ends.fasta -bed eccdna_ends
    awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
    awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
    paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq
    parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}
done

paste $(find . -maxdepth 1 -name "${INPUT_FILE}_shuffle_*_outputtable" | xargs -r ls -1 | tr "\n" " ") > ${INPUT_FILE}_shuffle_outputtable