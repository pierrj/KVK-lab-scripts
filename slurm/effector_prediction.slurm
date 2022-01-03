#!/bin/bash
#SBATCH --job-name=effector_prediction
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/effector_prediction

/global/scratch/users/pierrj/signalp-4.1/signalp -t euk -u 0.34 -U 0.34 -f short guy11_fungap_out_12_28_20_prot.faa > guy11_signalp_short


cat guy11_fungap_out_12_28_20_prot.faa | /global/scratch/users/pierrj/tmhmm-2.0c/bin/tmhmm > guy11_tmhmm_out