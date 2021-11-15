#!/bin/bash
#SBATCH --job-name=permutation_neverfound_genes_final
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

# bed file of gene locations
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
# scaffold and scaffold lengths
CHROM_SIZE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound

if [ -f "all.ecc_caller_out.regions.bed" ]; then
    rm all.ecc_caller_out.regions.bed
fi


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

## write all eccdna forming regions, minus lowq
while read sample; do
cd ${sample}
    cat ${sample}.ecc_caller_out.details.nolowq.txt >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound/all.ecc_caller_out.regions.bed
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound


## steps for rarefaction
## should be eccdnas not splitreads here technically
splitreadcount=$(wc -l  all.ecc_caller_out.regions.bed | awk '{print $1}')
splitreadcount_steps=$(wc -l  all.ecc_caller_out.regions.bed | awk '{print int($1/10)}')



for k in $(seq 0 99); do # 100 replicates
    if [ -f "genes_found_rarefaction_${k}" ]; then
        rm genes_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do ## go step by step for rarefaction
        shuf all.ecc_caller_out.regions.bed | head -n ${i} > subset.all.ecc_caller_out.regions.bed ## pick a random set of eccdnas
        ## count gene founds fully encompassed by eccdnas
        bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b subset.all.ecc_caller_out.regions.bed | awk -v OFS='\t' '{print $4, $5}' | awk '{ if ($2!=0) {print $1}}' | wc -l >> genes_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "genes_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > genes_found_rarefaction

for k in $(seq 0 99); do # 100 replicates
    if [ -f "random_found_rarefaction_${k}" ]; then
        rm random_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do ## go step by step for rarefaction
        # shuffle eccdnas across genome
        bedtools shuffle -i all.ecc_caller_out.regions.bed  -g ${CHROM_SIZE} | head -n ${i} > subset.all.ecc_caller_out.regions.bed
        ## count gene founds fully encompassed by randomly selected regions
        bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b subset.all.ecc_caller_out.regions.bed | awk -v OFS='\t' '{print $4, $5}' | awk '{ if ($2!=0) {print $1}}' | wc -l >> random_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "random_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > random_found_rarefaction