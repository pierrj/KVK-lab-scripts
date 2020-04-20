#!/bin/bash
#SBATCH --job-name=exonerate
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
exonerate --model protein2genome /global/scratch/users/pierrj/references/M.oryzae.EFCs.fasta /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no | grep "gene" | awk '$19>70' >  /global/scratch/users/pierrj/references/M.oryzae.EFCs.gff