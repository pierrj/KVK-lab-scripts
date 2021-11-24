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


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads

# bed file containing locations of te consensus sequences generated from wickersoft
TE_BEDFILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_nolengthfilter_individual_elements.bed
# mapfile with element names
ELEMENT_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots/transposon_mapfile

while read sample
do
    bamfile=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam
    # get split reads that could be attributed to ltr ecc circles
    /global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${sample} -b ${bamfile}
    split_read_file=candidate_ltr_srs.${sample}.bam
    # filter split reads to those attributed to ltr ecc circles
    /global/home/users/pierrj/git/bash/get_ltr_sr_circle_count_all_elements.sh -b ${TE_BEDFILE} \
        -m ${ELEMENT_MAPFILE} \
        -s ${split_read_file} -o ${sample}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

## normalize and average across bioreps and samples


# generate mapfile
/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

while read sample; do
    if [ -d "${sample}" ]; then
        rm -r ${sample}
    fi
    mkdir ${sample}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

if [ -f "results_mapfile" ]; then
    rm results_mapfile
fi

# names of results to grab
echo ltr_sr_cov_perfeature >> results_mapfile
echo ltr_ltr_sr_cov_perfeature >> results_mapfile
echo ltr_internal_sr_cov_perfeature >> results_mapfile
echo read_cov_perfeature >> results_mapfile
echo junction_sr_cov_perfeature >> results_mapfile

if [ -f "cov_perfeature_table" ]; then
    rm cov_perfeature_table
fi

awk '{print $1}' junction_sr_cov_perfeature_normalized > cov_perfeature_table

while read results; do
    # move results around
    while read sample; do
        cp ${sample}.${results} ${sample}/${sample}.${results}
    done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
    if [ -f "G3.normalization_table.${results}" ]; then
        rm G3.normalization_table.${results}
    fi
    # make normalization table
    while read sample; do
        ## normalization is a little tricky here but simplest is to just go by number of reads
        total=$( samtools view -c -F 4 -F 2048 /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam | awk '{print $1/1000000}')
        echo -e ${sample}'\t'${total} >> G3.normalization_table.${results}
    done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
    # normalize and average
    /global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
        -t G3_1A.${results} -n G3.normalization_table.${results} -y t
    /global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
    mv G3.normalized_binned ${results}_normalized
    paste cov_perfeature_table <(awk '{print $2}' ${results}_normalized) > tmptable
    mv tmptable cov_perfeature_table
done < results_mapfile


# make final table
cp cov_perfeature_table cov_feature_table

# these are all the same
# table for read normalization
cp G3.normalization_table.ltr_sr_cov_perfeature G3.normalization_table

## make cov_feature_table per sample, dont normalize though
while read sample; do
    if [ -f "${sample}.cov_perfeature_table" ]; then
        rm ${sample}.cov_perfeature_table
    fi
    awk '{print $1}' G3_1A/G3_1A.ltr_sr_cov_perfeature > ${sample}.cov_perfeature_table
    while read results; do
        paste ${sample}.cov_perfeature_table <( awk '{print $2}' ${sample}/${sample}.${results}) > ${sample}.tmptable
        mv ${sample}.tmptable ${sample}.cov_perfeature_table
    done < results_mapfile
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile