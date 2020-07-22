#!/bin/bash
#SBATCH --job-name=OsSR_noseqprep
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
mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample="OsSR_noseqprep"

cd ${SLURM_SUBMIT_DIR} ## submit from pipeline_tests
cd ${sample}

/global/home/users/pierrj/git/bash/generate_bam_file_genomecov_noseqprep.sh -g ${genome_bwa} \
 -1 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_R1.fastq -2 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_R2.fastq \
 -s ${sample} -t ${SLURM_NTASKS}

/global/home/users/pierrj/git/bash/generate_coverageplot.sh -m ${mapfile} -s ${sample} -c ${sample}.genomecoverage.bed

/global/home/users/pierrj/git/bash/call_ecc_regions.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b ${sample}.sorted.mergedandpe.bwamem.bam

/global/home/users/pierrj/git/bash/assign_confidence.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -c ${sample}.genomecoverage.bed