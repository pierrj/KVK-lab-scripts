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
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots

# list of chromosome names and sizes
genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
track=ltr_ecc_density
group=ltr_eccs

# genomic windows
bedtools makewindows -g ${CHROM_SIZES} -w 100 > guy11_genome_baoetal2017.100bpwindows

while read sample
do
# junction split reads corresponding to each window
bedtools coverage -a guy11_genome_baoetal2017.100bpwindows \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/ltr_eccs_splitreads -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.genomecoverage.bed
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

# table of read counts per sample
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

# make table of samples and numbers to use for normalization
/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.${group}.genomecoverage.bed -n ${group}.normalization_table -y t

# normalize and average coverage across all samples, technical replicates are averaged across biological replicates, biological replicates are then averaged
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

mv G3.normalized_binned G3.ltr_splitreads_coverage




MAPFILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
SAMPLE=G3
COVFILE=G3.ltr_splitreads_coverage

# filter to bins of interest (just a double check, to use with other organisms if necessary)
awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.normalized_binned.filtered
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
# rename bins
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.normalized_binned.filtered > ${SAMPLE}.normalized_binned.filtered.renamed



# name bins for input into manhattan plot script
awk -v OFS='\t' '{print $1, $2, $2+99, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > binned.${SAMPLE}.normalized_binned.filtered.renamed.bed
awk -v OFS='\t' '{print $1, $2, $3, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT''\t''SNP' > tmp.first_row
cat tmp.first_row ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table

cp /global/scratch/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_nolengthfilter.bed .

# make mapfile of transposon names
awk '{if ($1 !~ /^#/) {print substr(substr($4,16),0,length(substr($4,16))-1)}}' moryzae_high_conf_ltr_tes_nolengthfilter.bed | sort | uniq > transposon_mapfile

# make file of transposon locations for highlighting in manhattan plot
while read transposon; do
    awk -v t=${transposon} '$4 ~ t' moryzae_high_conf_ltr_tes_nolengthfilter.bed > ${transposon}.bed
    awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${transposon}.bed > renamed.${SAMPLE}.highlight_locs.bed
    bedtools intersect -wa -f 0.5 -a binned.${SAMPLE}.normalized_binned.filtered.renamed.bed -b renamed.${SAMPLE}.highlight_locs.bed | cut -f 4 > ${transposon}.snps.tohighlight
done < transposon_mapfile





### NOW WITH ALL SPLIT READS

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots_all_splitreads

# # list of chromosome names and sizes
genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
group=all_eccs

# genomic windows
bedtools makewindows -g ${CHROM_SIZES} -w 100 > guy11_genome_baoetal2017.100bpwindows

while read sample
do
# junction split reads corresponding to each window
bedtools coverage -a guy11_genome_baoetal2017.100bpwindows \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.splitreads.bed -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.genomecoverage.bed
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

# table of read counts per sample
while read sample; do
    total=$(wc -l /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.splitreads.bed | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${total} >> ${group}.normalization_table
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

while read sample; do
    if [ -d "${sample}" ]; then
        rm -r ${sample}
    fi
    mkdir ${sample}
    mv ${sample}.${group}.genomecoverage.bed ${sample}/${sample}.${group}.genomecoverage.bed
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

# make table of samples and numbers to use for normalization
/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots/sample_mapfile \
    -t G3_1A.${group}.genomecoverage.bed -n ${group}.normalization_table -y t

# normalize and average coverage across all samples, technical replicates are averaged across biological replicates, biological replicates are then averaged
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

mv G3.normalized_binned G3.all_splitreads_coverage


MAPFILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
SAMPLE=G3
COVFILE=G3.all_splitreads_coverage

# filter to bins of interest (just a double check, to use with other organisms if necessary)
awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.normalized_binned.filtered
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
# rename bins
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.normalized_binned.filtered > ${SAMPLE}.normalized_binned.filtered.renamed



# name bins for input into manhattan plot script
awk -v OFS='\t' '{print $1, $2, $2+99, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > binned.${SAMPLE}.normalized_binned.filtered.renamed.bed
awk -v OFS='\t' '{print $1, $2, $3, $1"_"$2}' ${SAMPLE}.normalized_binned.filtered.renamed > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames
echo -e 'CHROMOSOME''\t''BASE''\t''COUNT''\t''SNP' > tmp.first_row
cat tmp.first_row ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames > ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table

mv ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table ${SAMPLE}.normalized_binned.filtered.renamed.withbinnames.table.all_splitreads