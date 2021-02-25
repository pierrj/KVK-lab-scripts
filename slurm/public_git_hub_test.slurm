#!/bin/bash
#SBATCH --job-name=public_git_hub_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_caller_public_version_tests

cd ecc_caller
git pull
cd ..

export ECC_CALLER_PYTHON_SCRIPTS=ecc_caller/python_scripts/
export ECC_CALLER_PICARD=/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/

MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"

ecc_caller/generate_bam_file.sh -g ${GENOME_DB} -1 R1.fastq -2 R2.fastq -s output_name -t ${SLURM_NTASKS} -m ${MAPFILE}

ecc_caller/call_ecc_regions.sh -m ${MAPFILE} -s output_name -t ${SLURM_NTASKS} -b uniq.filtered.sorted.output_name.bam -q multimapped.filtered.name_sorted.output_name.bam

ecc_caller/assign_confidence.sh -m ${MAPFILE} -s output_name -t ${SLURM_NTASKS} -b no_secondary.filtered.sorted.output_name.bam -r output_name.confirmedsplitreads.bed