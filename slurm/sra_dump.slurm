#!/bin/bash
#SBATCH --job-name=sra_fastqdump
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/genome/2018_foster
# while read p; do
#   /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e 24 -O /global/scratch/users/pierrj/rnaseq/2019_luciano -t /global/scratch/users/pierrj/rnaseq/2019_luciano/tmp "$p"
# done < SRR_list.txt
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ERR2660591 -O /global/scratch/users/pierrj/genome/2018_foster
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e 24 -O /global/scratch/users/pierrj/genome/2018_foster -t /global/scratch/users/pierrj/genome/2018_foster/tmp /global/scratch/users/pierrj/genome/2018_foster/ERR2660591.sra
# mv ERR2660591_R1.fastq ERR2660591_R1.fastq.old
# mv ERR2660591_R2.fastq ERR2660591_R2.fastq.old
# awk '{ if (NR % 2 == 1) print $1, $2, $2".1", $3 ; else print $0}' ERR2660591_R1.fastq.old > ERR2660591_R1.fastq
# awk '{ if (NR % 2 == 1) print $1, $2, $2".1", $3 ; else print $0}' ERR2660591_R2.fastq.old > ERR2660591_R2.fastq