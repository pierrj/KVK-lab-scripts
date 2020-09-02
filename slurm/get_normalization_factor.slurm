#!/bin/bash
#SBATCH --job-name=get_normalization_factor
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
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A

slope_file="/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/slope_file"
spike_mapfile="/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/spike_mapfile"

if [ -f "${slope_file}" ]; then
    rm ${slope_file}
fi

/global/home/users/pierrj/git/bash/get_normalize_ratio.sh -m ${spike_mapfile} -s CT_1A -b CT_1A.sorted.mergedandpe.bwamem.bam -t ${SLURM_NTASKS}