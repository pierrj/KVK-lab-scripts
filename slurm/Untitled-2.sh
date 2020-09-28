#!/bin/bash
#SBATCH --job-name=bao_pfam_scan
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
cd /global/scratch/users/pierrj/pfam/bao_pfam_scan
source activate pfamscan

pfam_scan.pl -cpu 24 -outfile /global/scratch/users/pierrj/pfam/bao_pfam_scan/bao.pfamscan.txt -dir /global/scratch/users/pierrj/pfam -fasta /global/scratch/users/pierrj/references/Guy11_augustus_raw.aa

/global/scratch/users/pierrj/scripts/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p bao.pfamscan.txt  -e 0.001 -obao.pfamscan.kparse.txt
