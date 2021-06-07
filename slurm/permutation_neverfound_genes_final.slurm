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


GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
CHROM_SIZE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound

if [ -f "all.ecc_caller_out.regions.bed" ]; then
    rm all.ecc_caller_out.regions.bed
fi


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/


while read sample; do
cd ${sample}
    cat ${sample}.ecc_caller_out.details.nolowq.txt >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound/all.ecc_caller_out.regions.bed
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound

splitreadcount=$(wc -l  all.ecc_caller_out.regions.bed | awk '{print $1}')
splitreadcount_steps=$(wc -l  all.ecc_caller_out.regions.bed | awk '{print int($1/10)}')



for k in $(seq 0 99); do
    if [ -f "genes_found_rarefaction_${k}" ]; then
        rm genes_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do
        shuf all.ecc_caller_out.regions.bed | head -n ${i} > subset.all.ecc_caller_out.regions.bed
        bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b subset.all.ecc_caller_out.regions.bed | awk -v OFS='\t' '{print $4, $5}' | awk '{ if ($2!=0) {print $1}}' | wc -l >> genes_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "genes_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > genes_found_rarefaction

for k in $(seq 0 99); do
    if [ -f "random_found_rarefaction_${k}" ]; then
        rm random_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do
        bedtools shuffle -i all.ecc_caller_out.regions.bed  -g ${CHROM_SIZE} | head -n ${i} > subset.all.ecc_caller_out.regions.bed
        bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b subset.all.ecc_caller_out.regions.bed | awk -v OFS='\t' '{print $4, $5}' | awk '{ if ($2!=0) {print $1}}' | wc -l >> random_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "random_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > random_found_rarefaction