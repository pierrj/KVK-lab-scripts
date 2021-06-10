#!/bin/bash
#SBATCH --job-name=permutation_eccdna_forming_regions_details
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



if [ -f "all.ecc_caller_out.details.nolowq.txt" ]; then
    rm all.ecc_caller_out.details.nolowq.txt
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/


while read sample; do
cd ${sample}
    cat ${sample}.ecc_caller_out.details.nolowq.txt >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/all.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

splitreadcount=$(wc -l  all.ecc_caller_out.details.nolowq.txt | awk '{print $1}')
splitreadcount_steps=$(wc -l  all.ecc_caller_out.details.nolowq.txt | awk '{print int($1/10)}')



for k in $(seq 0 99); do
    if [ -f "regions_found_rarefaction_${k}" ]; then
        rm regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do
        shuf all.ecc_caller_out.details.nolowq.txt | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> details_regions_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "details_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > details_regions_found_rarefaction