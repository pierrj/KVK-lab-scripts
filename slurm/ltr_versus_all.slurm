#!/bin/bash
#SBATCH --job-name=ltr_versus_all
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons
rm ltr_versus_all
touch ltr_versus_all
while read sample;
do
    total=$(samtools view -c /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/guy11.sorted.allmapq.mergedandpe.${sample}_bwamem.bam)
    te_sum=$(awk '{if ($12 ~ /Gypsy/ || $9 ~ /Copia/) print $0}' ${sample}.teintersect.bed | awk '{sum+=$15} END {print sum}')
    echo "scale=4;($total-$nonefcs)/$total" | bc
    percent=$(echo "scale=4;$te_sum/$total" | bc)
    echo -e ${sample}'\t'$te_sum'\t'$total'\t'$percent >> ltr_versus_all
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
while read sample;
do
    total=$(samtools view -c /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/${sample}/guy11.sorted.allmapq.mergedandpe.${sample}_bwamem.bam)
    te_sum=$(awk '{if ($12 ~ /Gypsy/ || $9 ~ /Copia/) print $0}' ${sample}.teintersect.bed | awk '{sum+=$15} END {print sum}')
    echo "scale=4;($total-$nonefcs)/$total" | bc
    percent=$(echo "scale=4;$te_sum/$total" | bc)
    echo -e ${sample}'\t'$te_sum'\t'$total'\t'$percent >> ltr_versus_all
done < /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/mapfile