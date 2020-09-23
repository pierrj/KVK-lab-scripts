#!/bin/bash
#SBATCH --job-name=split_reads_per_ltr_per_sample_normalization
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

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m mapfile_IF

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t IF_1A.sr_count_per_element -n SRs.no_orientation.allmapq.IF_1A.bam -y e

cd results/ltr_splitreads_readnormalized

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ../../mapfile_for_normalize_and_average_ltrs_totalreads -f 1000000 -b 1 -c 2 -n f