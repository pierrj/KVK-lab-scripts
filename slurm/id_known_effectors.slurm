#!/bin/bash
#SBATCH --job-name=id_known_effectors
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/diamond/guy11_avr_locs

# proteome from fungap
../diamond makedb --in guy11_fungap_out_12_28_20_prot.faa -d guy11_fungap_12_28_20_prot

# avrs from Kim et al 2019
../diamond blastp -d guy11_fungap_12_28_20_prot -q M.oryzae.AVRs.fasta -o matches.tsv


## filter out identity less than 50%, then make sure each fungap gene only has the best match assigned to it
awk '$3>=50.0' matches.tsv | awk '$2 != "gene_09984.t1" && $1 != "AvrPita2" && $1 != "Pwl3" && !($1 == "AvrPita3" && $2 == "gene_10084.t1")' \
 | awk -v OFS='\t' '{print $1, $2}' > known_effectors_fungap_names

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common

cp /global/scratch/users/pierrj/diamond/guy11_avr_locs/known_effectors_fungap_names .