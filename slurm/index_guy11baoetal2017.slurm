#!/bin/bash
#SBATCH --job-name=mappbsim_restrictive
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --time=23:59:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/analysis/pbsim_mapped_restrictive
source activate python=2.7
pbmm2 index -k 28 ../../reference_genomes/guy11_genome_baoetal2017.fasta guy11_genomebaoetal2017_restrictive.mmi
pbmm2 align --preset CCS --sort --log-file pbsimreads_mapped_restrictive.log --log-level INFO guy11_genomebaoetal2017_restrictive.mmi ../../rawdata/pbsim/pooled_pbsimreads.fastq pbsimreads_mapped_restrictive.bam
source deactivate
