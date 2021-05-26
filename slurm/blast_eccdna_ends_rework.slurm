#!/bin/bash
#SBATCH --job-name=blast_eccDNA_ends_rework
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/homology/

INPUT_FILE=shuffle_eccs.bed
OUTPUT_TABLE=random_homology_90_perc
GENOME=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

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

## big issue here with skipped getfastas!
awk '$1 !~ /^>/' eccdna_starts.fasta > nonames_eccdna_starts.seq
awk '$1 !~ /^>/' eccdna_ends.fasta > nonames_eccdna_ends.seq
paste nonames_eccdna_starts.seq nonames_eccdna_ends.seq > eccdnas_starts_and_ends.seq

parallel --colsep '\t' 'blastn -query <(echo -e \>test"\n"{1}) -subject <(echo -e \>test"\n"{2}) -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 90 -outfmt  "6 length" | awk "{print}; END {if (NR<1) {print 0}}"' < eccdnas_starts_and_ends.seq > ${OUTPUT_TABLE}