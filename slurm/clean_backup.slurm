#!/bin/bash
#SBATCH --job-name=clean_backup
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/projects/vector_kvklab/pierrj/eccdna_backup/eccDNA

# find . -name "*.fastq" -type f -delete

# find . -name "*.sra" -type f -delete

# find . -name "*.bam" -type f -delete

# find . -name "*.sam" -type f -delete

rm -r rerunning_stuff_final/pipeline_qc/rarefaction_tech_reps/

rm -r rerunning_stuff_final/pipeline_qc/pacbio/all_false_positives

cd /global/scratch/users/pierrj/space_snif

du -h /global/scratch/projects/vector_kvklab/pierrj/eccdna_backup/eccDNA > eccdna_dir_space_sniff_post_rm_bam_files_post_rm_techreps_fps