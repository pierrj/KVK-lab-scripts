#!/bin/bash
#SBATCH --job-name=rice_pipelinetest
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/anygenome_rice_test
sample="RC_1A"
genome="/global/scratch/users/pierrj/references/guy11bao_7015mito_spikeins_tweaked_osativa_worganelles_bwa"
read_path="/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/"
mapfile="/global/scratch/users/pierrj/references/ORSA_w_organelles.contignames"

#/global/home/users/pierrj/git/bash/generate_bam_file_genomecov.sh -g ${genome} \
#  -1 ${read_path}${sample}/${sample}_R1.fastq -2 ${read_path}${sample}/${sample}_R2.fastq \
#  -s ${sample} -t ${SLURM_NTASKS}

#/global/home/users/pierrj/git/bash/generate_coverageplot.sh -m ${mapfile} -s ${sample} -c ${sample}.genomecoverage.bed

echo "yo"

/global/home/users/pierrj/git/bash/call_ecc_regions.sh -m ${mapfile} \
    -s ${sample} -t 10 -b ${sample}.sorted.mergedandpe.bwamem.bam

echo "whatsup"

/global/home/users/pierrj/git/bash/asign_confidence.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -c ${sample}.genomecoverage.bed