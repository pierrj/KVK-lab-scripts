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

sample=G3_2
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=G3_3
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=CT_2
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=CT_3
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=IF_2
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=IF_3
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=PQ_2
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done

sample=PQ_3
rm simulations_${sample}
touch simulations_${sample}
for i in $(seq 0 1 1000); do
    bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
    nonefcs=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.gff -b shuffled.${sample}.bed | awk '$5!=0' | wc -l)
    echo "scale=4;($total-$nonefcs)/$total" | bc >> simulations_${sample}
done