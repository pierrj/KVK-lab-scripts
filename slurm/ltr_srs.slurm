#!/bin/bash
#SBATCH --job-name=ltr_srs
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
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A
mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
sample="RC_1A"
sorted_bamfile="RC_1A.sorted.mergedandpe.bwamem.bam"

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -m ${mapfile} -s ${sample} -b ${sorted_bamfile}