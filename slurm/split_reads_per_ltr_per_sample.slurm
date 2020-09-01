#!/bin/bash
#SBATCH --job-name=split_reads_per_ltr_per_sample
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
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/

while read sample
do
cd ${sample}
bamfile=filtered.sorted.${sample}.bam
/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${sample} -b ${bamfile}
split_read_file=SRs.no_orientation.allmapq.${sample}.bam
/global/home/users/pierrj/git/bash/get_ltr_sr_circle_count_all_elements.sh -b /global/scratch/users/pierrj/repeatmasker/moryzae_final/copia_int_renamed.only_LTRs_of_interest.identitycutoff.bed \
    -m /global/scratch/users/pierrj/repeatmasker/moryzae_final/mapfile_elements \
    -s ${split_read_file} -o ${sample}
cd ..
done < mapfile_CT_PQ