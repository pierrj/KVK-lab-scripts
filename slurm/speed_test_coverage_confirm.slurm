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

mapfile="/local/job${SLURM_JOBID}/ORSA_no_organelles.contignames"
sample="RC_1A"
filtered_bamfile="/local/job${SLURM_JOBID}/filtered.sorted.RC_1A.bam"
bamfile="/local/job${SLURM_JOBID}/RC_1A.sorted.mergedandpe.bwamem.bam"


cd /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_no_organelles.contignames .
cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/filtered.sorted.RC_1A.bam .
cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A.sorted.mergedandpe.bwamem.bam .

mkdir nodb

cd nodb

cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/speed_test_coverage_confirm/nodb/parallel.confirmed .

source activate bedtools

date +"%T"

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b ${filtered_bamfile}

date +"%T"

source deactivate

cd ..

mkdir yesdb

cd yesdb

cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/speed_test_coverage_confirm/nodb/parallel.confirmed .

date +"%T"

bedtools genomecov -d -ibam ${bamfile} > ${sample}.genomecoverage.bed

/global/home/users/pierrj/git/bash/assign_confidence.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -c ${sample}.genomecoverage.bed

date +"%T"