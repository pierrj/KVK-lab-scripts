#!/bin/bash
#SBATCH --job-name=permute_test_deletions
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/deletion_mummer_plots/guy11_eccdnas

for i in {1..1000}; do 
    bedtools shuffle -i sorted.all.confirmedsplitreads.bed -g guy11_genome_baoetal2017.chromsizes > shuffled
    python /global/home/users/pierrj/git/python/check_hotspot_ecc_and_deletions.py uniq_deletions.bed shuffled
done

for i in {1..1000}; do 
    bedtools shuffle -i uniq_deletions.bed -g guy11_genome_baoetal2017.chromsizes > shuffled
    python /global/home/users/pierrj/git/python/check_overlap_ecc_and_deletions.py shuffled sorted.uniq.all.confirmedsplitreads.bed
done

for i in {1..1000}; do 
    bedtools shuffle -i uniq_deletions_real.bed -g guy11_genome_baoetal2017.chromsizes > shuffled
    python /global/home/users/pierrj/git/python/check_hotspot_ecc_and_deletions.py shuffled sorted.all.confirmedsplitreads.bed
done


for i in {1..1000}; do 
    bedtools shuffle -i uniq_small_deletions.bed -g guy11_genome_baoetal2017.chromsizes > shuffled
    python /global/home/users/pierrj/git/python/check_overlap_ecc_and_deletions.py shuffled sorted.uniq.all.confirmedsplitreads.bed
done