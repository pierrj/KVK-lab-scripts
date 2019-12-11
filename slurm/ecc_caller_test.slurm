#!/bin/bash
#SBATCH --job-name=ecc_caller_test
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/G3_1A_tests
ipcluster start -n $SLURM_NTASKS --cluster-id="slurm-${SLURM_JOBID}" &
sleep 45
ipython /global/home/users/pierrj/git/python/ecc_caller.py samechromosome.exactlytwice.reverseread1.G3_1A_bwamem.bam discordantmappedreads.oppositefacing.bed genomecoverage.mergedandpe.G3_1A_bwamem.bed G3_1A
ipcluster stop --cluster-id="slurm-${SLURM_JOBID}"