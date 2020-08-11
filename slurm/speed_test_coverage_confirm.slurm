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

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
sample="RC_1A"
filtered_bamfile="/global/scratch/users/pierrj/eccDNA/pipeline_tests/RC_1A_noseqprep/filtered.sorted.RC_1A_noseqprep.bam"
bamfile="/global/scratch/users/pierrj/eccDNA/pipeline_tests/RC_1A_noseqprep/RC_1A_noseqprep.sorted.mergedandpe.bwamem.bam"

cd nodb 

source activate bedtools

date +"%T"

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b ${filtered_bamfile}

date +"%T"

# source deactivate

# cd ..

# cd yesdb

# date +"%T"

# bedtools genomecov -d -ibam ${bamfile} > ${sample}.genomecoverage.bed

# /global/home/users/pierrj/git/bash/assign_confidence.sh -m ${mapfile} \
#     -s ${sample} -t ${SLURM_NTASKS} -c ${sample}.genomecoverage.bed

# date +"%T"