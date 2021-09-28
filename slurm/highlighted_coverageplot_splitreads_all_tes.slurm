#!/bin/bash
#SBATCH --job-name=highlighted_coverageplot_splitreads_all_tes
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
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots

genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
track=ltr_ecc_density
group=ltr_eccs

bedtools makewindows -g ${CHROM_SIZES} -w 100 > guy11_genome_baoetal2017.100bpwindows

while read sample
do
bedtools coverage -a guy11_genome_baoetal2017.100bpwindows \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/ltr_eccs_splitreads -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.genomecoverage.bed
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

while read sample; do
    total=$(wc -l /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/ltr_eccs_splitreads | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${total} >> ${group}.normalization_table
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

while read sample; do
    if [ -d "${sample}" ]; then
        rm -r ${sample}
    fi
    mkdir ${sample}
    mv ${sample}.${group}.genomecoverage.bed ${sample}/${sample}.${group}.genomecoverage.bed
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.${group}.genomecoverage.bed -n ${group}.normalization_table -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

mv G3.normalized_binned G3.ltr_splitreads_coverage




MAPFILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
SAMPLE=G3
COVFILE=G3.ltr_splitreads_coverage

awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.normalized_binned.filtered
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.normalized_binned.filtered > ${SAMPLE}.normalized_binned.filtered.renamed



## WHAT DOES THIS DO THOUGH MAYBE HAS TO DO WITH BEDTOOLS INTERSECT
awk -v OFS='\t' '{print $1, $2, $2+99, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > binned.${SAMPLE}.normalized_binned.filtered.renamed.bed
awk -v OFS='\t' '{print $1, $2, $3, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT''\t''SNP' > tmp.first_row
cat tmp.first_row ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table

cp /global/scratch/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_nolengthfilter.bed .

awk '{if ($1 !~ /^#/) {print substr(substr($4,16),0,length(substr($4,16))-1)}}' moryzae_high_conf_ltr_tes_nolengthfilter.bed | sort | uniq > transposon_mapfile

while read transposon; do
    awk -v t=${transposon} '$4 ~ t' moryzae_high_conf_ltr_tes_nolengthfilter.bed > ${transposon}.bed
    awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${transposon}.bed > renamed.${SAMPLE}.highlight_locs.bed
    bedtools intersect -wa -f 0.5 -a binned.${SAMPLE}.normalized_binned.filtered.renamed.bed -b renamed.${SAMPLE}.highlight_locs.bed | cut -f 4 > ${transposon}.snps.tohighlight
done < transposon_mapfile