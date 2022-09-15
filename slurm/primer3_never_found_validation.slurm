#!/bin/bash
#SBATCH --job-name=primer3_never_found_validation
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

NUM_PRIMER_PAIRS=2
SLOP=250

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/neverfound_pcr/

if [ -f "absent_gene_locs.bed" ]; then
    rm absent_gene_locs.bed
fi

echo -e MQOP01000004.1"\t"2519225"\t"2523019"\t"MYO1 >> absent_gene_locs.bed
echo -e MQOP01000005.1"\t"1955859"\t"1957619"\t"actin >> absent_gene_locs.bed

bedtools slop -i absent_gene_locs.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes -b $SLOP > absent_gene_locs_slop.bed

bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -bed absent_gene_locs_slop.bed > absent_gene_locs_slop.fasta

while read -r ONE; do read -r TWO; echo "SEQUENCE_TEMPLATE=${TWO}"; done < absent_gene_locs_slop.fasta > absent_gene_locs_slop.templates

awk -v OFS='\t' '{print "SEQUENCE_ID="$4}' absent_gene_locs_slop.bed > absent_gene_locs_slop.names

a=($(wc -l absent_gene_locs_slop.names))


if [ -f "tmp.primer_params" ]; then
    rm tmp.primer_params
fi

# awk '{print $3-$2}' absent_gene_locs_slop.bed

echo "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,200,4094,200" >> tmp.primer_params
echo "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,200,2060,200" >> tmp.primer_params

if [ -f "tmp.primer_num" ]; then
    rm tmp.primer_num
fi

echo "PRIMER_NUM_RETURN=$NUM_PRIMER_PAIRS" >> tmp.primer_num
echo "PRIMER_NUM_RETURN=$NUM_PRIMER_PAIRS" >> tmp.primer_num

if [ -f "tmp.product_size" ]; then
    rm tmp.product_size
fi

echo "PRIMER_PRODUCT_SIZE_RANGE=100-6000" >> tmp.product_size
echo "PRIMER_PRODUCT_SIZE_RANGE=100-6000" >> tmp.product_size

if [ -f "tmp.equal_signs" ]; then
    rm tmp.equal_signs
fi

echo "=" >> tmp.equal_signs
echo "=" >> tmp.equal_signs

paste -d"\n" absent_gene_locs_slop.names absent_gene_locs_slop.templates tmp.primer_params tmp.primer_num tmp.product_size tmp.equal_signs > absent_gene_locs_slop.boulder

/global/scratch/users/pierrj/scripts/primer3/src/primer3_core \
    --p3_settings_file=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/neverfound_pcr/PMJ_settings_neverfound_primers \
    absent_gene_locs_slop.boulder > tmp.primer3.boulder

grep -oP '(PRIMER_.*._.*._SEQUENCE=)\K.*' tmp.primer3.boulder > tmp.primer_seqs

awk '{print $4}' absent_gene_locs.bed > tmp.gene_names_awk

if [ -f "tmp.primer_names" ]; then
    rm tmp.primer_names
fi

# grab primer names
while read gene; do
    for i in $(seq 1 $NUM_PRIMER_PAIRS); do
        echo ${gene}_F${i} >> tmp.primer_names
        echo ${gene}_R${i} >> tmp.primer_names
    done
done < tmp.gene_names_awk

if [ -f "tmp.primer_names_product_size" ]; then
    rm tmp.primer_names_product_size
fi

# grab product sizes
while read gene; do
    for i in $(seq 1 $NUM_PRIMER_PAIRS); do
        echo ${gene}_F${i}_${gene}_R${i} >> tmp.primer_names_product_size
    done
done < tmp.gene_names_awk

## make output file for export into IDT for ordering
paste -d';' tmp.primer_names tmp.primer_seqs > neverfound.idt

grep -oP '(PRIMER_PAIR_.*._PRODUCT_SIZE=)\K.*' tmp.primer3.boulder | paste tmp.primer_names_product_size - > neverfound.product_sizes