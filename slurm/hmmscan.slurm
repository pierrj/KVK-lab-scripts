#!/bin/bash
#SBATCH --job-name=pfam_scan
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
cd /global/scratch/users/pierrj/pfam/G3_BRAKER_final
module load hmmer

hmmscan --cpu 24 --domtblout G3.hmmscan.txt ../Pfam-A.hmm /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/genes_on_eccs/G3_all.common.genes.fasta

hmmscan --cpu 24 --domtblout WG.hmmscan.txt ../Pfam-A.hmm /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.CDS.fasta