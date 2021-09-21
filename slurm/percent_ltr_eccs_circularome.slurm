#!/bin/bash
#SBATCH --job-name=percent_ltr_eccs_circularome
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

if [ -f "ltr_eccs_count" ]; then
    rm ltr_eccs_count
fi

if [ -f "ltr_eccs_splitreads_count" ]; then
    rm ltr_eccs_splitreads_count
fi

if [ -f "ecc_region_count" ]; then
    rm ecc_region_count
fi

if [ -f "junction_sr_count" ]; then
    rm junction_sr_count
fi

while read sample; do
cd ${sample}
    wc -l ltr_eccs | awk '{print $1}' >> ../ltr_eccs_count
    wc -l ltr_eccs_splitreads | awk '{print $1}' >> ../ltr_eccs_splitreads_count
    wc -l ${sample}.ecc_caller_out.details.nolowq.txt | awk '{print $1}' >> ../ecc_region_count
    wc -l ${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}' >> ../junction_sr_count
cd ..
done < mapfile

awk '{sum+=$1} END {print sum}' ltr_eccs_count
awk '{sum+=$1} END {print sum}' ltr_eccs_splitreads_count
awk '{sum+=$1} END {print sum}' ecc_region_count
awk '{sum+=$1} END {print sum}' junction_sr_count
