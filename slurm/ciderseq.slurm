#!/bin/bash
#SBATCH --job-name=ciderseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all/G3_1A

module remove python

source activate ciderseq

module load muscle

SAMPLE=G3_1A

THREADS=${SLURM_NTASKS}

python3 /global/scratch/users/pierrj/ciderseq2/cs-tools.py split --mode deconcat --format fastq --numjobs ${THREADS} --cluster "srun --ntasks=1" /global/scratch/users/pierrj/ciderseq2/ciderseq_config.json ${SAMPLE}_pacbio.fastq
python3 /global/scratch/users/pierrj/ciderseq2/cs-tools.py join --clean /global/scratch/users/pierrj/ciderseq2/ciderseq_config.json ${SAMPLE}_pacbio.fastq
python3 /global/scratch/users/pierrj/ciderseq2/cs-tools.py plot /global/scratch/users/pierrj/ciderseq2/ciderseq_config.json ${SAMPLE}_pacbio.fastq