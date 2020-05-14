#!/bin/bash
#SBATCH --job-name=pfam_scan_WG
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
cd /global/scratch/users/pierrj/pfam/
module load hmmer
hmmscan --cpu 24 --domtblout WG.hmmscan.txt Pfam-A.hmm /global/scratch/users/pierrj/references/Guy11_augustus_raw.aa

# rm /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/geneids.G3_1and2and3.fasta
# touch /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/geneids.G3_1and2and3.fasta
# while read geneid;
# do
#     grep -A1 ${geneid}.t1 Guy11_augustus_raw.aa >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/geneids.G3_1and2and3.fasta
# done < /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/geneids.G3_1and2and3
