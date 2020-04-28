#!/bin/bash
#SBATCH --job-name=efc_enrichment_redistributed
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment
rm redistributed_enrichment
touch redistributed_enrichment
while read sample;
do
    nonefcs=$(bedtools intersect -f 1 -wa -c -a moryzae_all.redistributed.gff -b ../merging/parallel.confirmed.${sample}.bed | awk '$5!=0' | awk '$4 ~ /^ID/' | wc -l)
    total=$(bedtools intersect -f 1 -wa -c -a moryzae_all.redistributed.gff -b ../merging/parallel.confirmed.${sample}.bed | awk '$5!=0' | wc -l)
    echo -e ${sample}"\t"$(echo "scale=4;($total-$nonefcs)/$total" | bc) >> redistributed_enrichment
done < mapfile