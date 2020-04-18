#!/bin/bash
#SBATCH --job-name=demux_ccs
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/SMRTcell1
source activate python=2.7
ccs --reportFile=ccs_report_smrtcell1_singlenode.txt subreads.bam ccs.bam
lima --different --dump-removed --split-bam-named --ccs ccs.bam /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/barcodes.fasta demux.bam