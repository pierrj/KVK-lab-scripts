#!/bin/bash
#SBATCH --job-name=simulated_eccDNAs_all
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/distribution

# while read sample;
# do
#     rm simulations_${sample}
#     touch simulations_${sample}
#     awk -v OFS='\t' '{print $1, $2, $3}' ../merging/parallel.confirmed.${sample}.bed | sort | uniq  > parallel.confirmed.${sample}.unique.bed
#     for i in $(seq 0 1 1000); do
#         bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
#         nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
#         total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
#         echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
#     done
# done < mapfile


rm observed_expected
touch observed_expected

while read sample;
do
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b parallel.confirmed.${sample}.unique.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b parallel.confirmed.${sample}.unique.bed | awk '$5!=0' | wc -l)
    observed=$(echo "scale=4;($total-$nonefcs)/$total" | bc)
    expected=$(awk '{sum+=$1} END {print sum/NR}' simulations_${sample})
    echo -e ${sample}'\t'$observed'\t'$expected >> observed_expected
done < mapfile