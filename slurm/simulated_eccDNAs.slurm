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
rm observed_expected
touch observed_expected
while read sample;
do
    rm simulations_${sample}
    touch simulations_${sample}
    awk -v OFS='\t' '{print $1, $2, $3}' ../merging/parallel.confirmed.${sample}.bed | sort | uniq  > parallel.confirmed.${sample}.unique.bed
    efc_containing=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_efcsandavrs_trimmed.uniquefiltered.gff -b parallel.confirmed.${sample}.unique.bed | awk '{sum+=$5;} END{print sum;}')
    total=$(wc -l parallel.confirmed.${sample}.unique.bed | awk '{print $1}')
    observed_proportion=$(echo "scale=4;$efc_containing/$total" | bc)
    for i in $(seq 0 1 1000); do
        bedtools shuffle -i parallel.confirmed.${sample}.unique.bed -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes > shuffled.${sample}.bed
        eccs_containing=$(bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_efcsandavrs_trimmed.uniquefiltered.gff -b shuffled.${sample}.bed | awk '{sum+=$5;} END{print sum;}')
        eccs_total=$(wc -l shuffled.${sample}.bed | awk '{print $1}')
        echo "scale=4;($eccs_containing)/$eccs_total" | bc >> simulations_${sample}
    done
    expected_proportion=$(awk '{sum+=$1} END {print sum/NR}' simulations_${sample})
    echo -e ${sample}'\t'$efc_containing'\t'$total'\t'$observed_proportion'\t'$expected_proportion >> observed_expected
done < mapfile