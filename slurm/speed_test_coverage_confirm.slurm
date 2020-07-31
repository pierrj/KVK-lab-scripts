#!/bin/bash
#SBATCH --job-name=speed_test_coverage_confirm
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/speed_test_coverage_confirm

mapfile="/global/scratch/users/pierrj/references/ORSA_w_organelles.contignames"
sample="RC_1A"
filtered_bamfile="/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/filtered.sorted.RC_1A.bam"
bamfile="/global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A.sorted.mergedandpe.bwamem.bam"

cd nodb

date +"%T"

source activate bedtools

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b ${filtered_bamfile}

source deactivate

cd ../yesdb

date +"%T"

bedtools genomecov -d -ibam ${bamfile} > ${SAMPLE}.genomecoverage.bed

/global/home/users/pierrj/git/bash/assign_confidence.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -c ${sample}.genomecoverage.bed

date +"%T"